use std::{collections::HashMap, str};

use anyhow::{bail, Context, Result};
use rust_htslib::bcf::{record::GenotypeAllele, Header, Record, Writer};

use crate::repeat::TandemRepeat;

/// Supported filter tags. If more variants are added,
/// they should also be added to `VCF_FILTER_LINES` in [`crate::io::vcf`]
#[derive(Debug)]
pub enum VcfFilter {
    Pass,
    Undef, // acts as a flag to skip locus without reporting reason in final VCF file
    DpZero,
    DpOor,
    CnZero,
    CnOor,
    CnMissing,
    AmbGt,
}

impl VcfFilter {
    pub fn name(&self) -> &str {
        match self {
            VcfFilter::Pass => "PASS",
            VcfFilter::Undef => "UNDEF",
            VcfFilter::DpZero => "DPZERO",
            VcfFilter::DpOor => "DPOOR",
            VcfFilter::CnZero => "CNZERO",
            VcfFilter::CnOor => "CNOOR",
            VcfFilter::CnMissing => "CNMISSING",
            VcfFilter::AmbGt => "AMBGT",
        }
    }
}

pub fn get_info_int(record: &Record, tag: &str) -> Result<i64> {
    let res = record
        .info(tag.as_bytes())
        .integer()
        .with_context(|| format!("Failed to parse info field '{tag}'"))?;

    if let Some(res) = res {
        Ok(res[0] as i64)
    } else {
        bail!("Info field '{tag}' was empty")
    }
}

pub fn get_info_str(record: &Record, tag: &str) -> Result<String> {
    let res = record
        .info(tag.as_bytes())
        .string()
        .with_context(|| format!("Failed to parse info field '{tag}'"))?;

    if let Some(res) = res {
        let s = str::from_utf8(res[0]).context("Error parsing VCF string field")?;
        Ok(s.to_string())
    } else {
        bail!("Info field '{tag}' was empty")
    }
}

#[allow(dead_code)]
fn get_format_int(record: &Record, tag: &str, sample_idx: usize) -> Result<Option<usize>> {
    let res = record.format(tag.as_bytes()).integer();
    if let Ok(res) = res {
        let sample_val = res[sample_idx];

        Ok(Some(sample_val[0] as usize))
    } else {
        Ok(None)
    }
}

fn get_format_str(record: &Record, tag: &str, sample_idx: usize) -> Result<Option<String>> {
    let res = record.format(tag.as_bytes()).string();

    if let Ok(res) = res {
        let sample_val = res[sample_idx];
        let s = str::from_utf8(sample_val).context("Error parsing VCF string field")?;

        Ok(Some(s.to_string()))
    } else {
        Ok(None)
    }
}

pub fn allele_lens_from_record(
    record: &Record,
    sample_idx: usize,
) -> Result<Option<HashMap<i64, f32>>> {
    let allele_freqs = get_format_str(record, "FREQS", sample_idx)?;
    if let Some(allele_freqs) = allele_freqs {
        let allele_freqs: Vec<&str> = allele_freqs.split('|').collect();

        let mut freq_map: HashMap<i64, f32> = HashMap::new();
        for i in allele_freqs.iter() {
            let split: Vec<&str> = i.split(',').collect();
            if split.len() != 2 {
                bail!("Encountered improperly formatted FREQS format field in VCF");
            }
            let allele_len = split[0].parse::<i64>()?;
            let allele_freq = split[1].parse::<f32>()?;
            freq_map.insert(allele_len, allele_freq);
        }

        Ok(Some(freq_map))
    } else {
        Ok(None)
    }
}

pub fn repeat_to_bcf_record(repeat: &TandemRepeat, writer: &Writer) -> Result<Record> {
    // Create record for repeat and add reference information
    let mut record = writer.empty_record();
    add_info_fields(&mut record, &writer, repeat)?;

    // Using the estimated genotype for this locus, create the alleles and genotype string to add to the VCF
    add_alleles_genotypes(&mut record, repeat)?;

    // Add additional information that was extracted from the alignment to the VCF af FORMAT fields
    add_format_fields(&mut record, repeat)?;

    Ok(record)
}

/// Add information about how the `tr_region` looks in the reference genome to the
/// VCF `record`. We need the `vcf` struct to get the contig id from the VCF header.
fn add_info_fields(record: &mut Record, vcf: &Writer, repeat: &TandemRepeat) -> Result<()> {
    let context = || {
        format!(
            "Error setting INFO field value for {}",
            repeat.reference_info.get_fetch_definition_s()
        )
    };

    let rid = vcf
        .header()
        .name2rid(repeat.reference_info.seqname.as_bytes())
        .with_context(context)?;
    let ref_len = repeat.reference_info.get_reference_len() as usize;
    record.set_rid(Some(rid));
    record.set_pos(repeat.reference_info.start);

    record
        .push_info_integer(b"END", &[repeat.reference_info.end as i32])
        .with_context(context)?;
    record
        .push_info_string(b"RU", &[repeat.reference_info.unit.as_bytes()])
        .with_context(context)?;
    record
        .push_info_integer(b"PERIOD", &[repeat.reference_info.period as i32])
        .with_context(context)?;
    record
        .push_info_integer(b"REF", &[ref_len as i32])
        .with_context(context)?;

    Ok(())
}

/// Add genotyping results stored on `tr_region` to the VCF `record`.
fn add_alleles_genotypes(record: &mut Record, repeat: &TandemRepeat) -> Result<()> {
    let ref_len = repeat.reference_info.get_reference_len() as usize;
    let ref_allele = repeat.reference_info.unit.repeat(ref_len);
    let mut alleles = vec![ref_allele];
    let mut genotype: Vec<GenotypeAllele> = Vec::new();

    match &repeat.genotype {
        None => genotype.push(GenotypeAllele::UnphasedMissing), // No genotype estimated for repeat
        Some(gts) => {
            // gts has form &(allele_length: i64, n_times_estimated: f32)
            for gt in gts {
                let len = gt.0 as i32;
                let n_estimated = gt.1 as i32;
                if len == ref_len as i32 {
                    // Reference allele always has to be first
                    // Add this genotype as many times as it was estimated to be present
                    for _ in 0..n_estimated {
                        genotype.insert(0, GenotypeAllele::Unphased(0));
                    }
                    continue;
                }
                // Add this genotype as many times as it was estimated to be present
                for _ in 0..n_estimated {
                    genotype.push(GenotypeAllele::Unphased((alleles.len()) as i32));
                }
                let allele = repeat.reference_info.unit.repeat(len as usize);
                alleles.push(allele);
            }
        }
    }

    let alleles: Vec<&[u8]> = alleles.iter().map(std::string::String::as_bytes).collect();
    record.set_alleles(&alleles).with_context(|| {
        format!(
            "Error setting alleles for {}",
            repeat.reference_info.get_fetch_definition_s()
        )
    })?;
    record.push_genotypes(&genotype).with_context(|| {
        format!(
            "Error setting genotype for {}",
            repeat.reference_info.get_fetch_definition_s()
        )
    })?;
    Ok(())
}

/// Add additional information stored on `tr_region` to the VCF `record` format fields
fn add_format_fields(record: &mut Record, repeat: &TandemRepeat) -> Result<()> {
    let context = || {
        format!(
            "Error setting FORMAT field value for {}",
            repeat.reference_info.get_fetch_definition_s()
        )
    };

    record
        .push_format_string(b"FT", &[repeat.filter.name().as_bytes()])
        .with_context(context)?;

    if !matches!(repeat.filter, VcfFilter::CnMissing) {
        record
            .push_format_integer(b"CN", &[repeat.copy_number as i32])
            .with_context(context)?;
    }

    if let Some(depth) = repeat.get_n_mapped_reads() {
        record
            .push_format_integer(b"DP", &[depth as i32])
            .with_context(context)?;
    } else {
        record
            .push_format_integer(b"DP", &[0])
            .with_context(context)?;
    }

    let mut allele_freqs = repeat.allele_freqs_as_tuples();
    allele_freqs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let allele_freqs: Vec<String> = allele_freqs
        .iter()
        .map(|i| format!("{},{}", i.0, i.1))
        .collect();
    let allele_freqs = allele_freqs.join("|");
    record
        .push_format_string(b"FREQS", &[allele_freqs.as_bytes()])
        .with_context(context)?;

    let gt_as_allele_lens: Vec<String> = repeat
        .gt_as_allele_lengths()
        .iter()
        .map(std::string::ToString::to_string)
        .collect();
    let gt_as_allele_lens = gt_as_allele_lens.join(",");
    record
        .push_format_string(b"REPLEN", &[gt_as_allele_lens.as_bytes()])
        .with_context(context)?;

    Ok(())
}

/// Construct VCF a header. First, include information about the target contigs, followed by the [`VCF_INFO_LINES`], [`VCF_FORMAT_LINES`].
/// Then, write TR variant calls to stdout.
pub fn make_vcf_header(targets: &[String], lengths: &[u64], sample_name: &str) -> Header {
    let mut header = Header::new();

    for (target, length) in targets.iter().zip(lengths.iter()) {
        let header_contig_line = format!(r#"##contig=<ID={target},length={length}>"#);
        header.push_record(header_contig_line.as_bytes());
    }
    for header_info_line in VCF_INFO_LINES {
        header.push_record(header_info_line);
    }
    for header_format_line in VCF_FORMAT_LINES {
        header.push_record(header_format_line);
    }
    for header_filter_line in VCF_FILTER_LINES {
        header.push_record(header_filter_line);
    }
    header.push_sample(sample_name.as_bytes());

    header
}

/// The VCF info lines to be included in the header. See [`make_vcf_header`].
const VCF_INFO_LINES: &[&[u8]] = &[
    br#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of reference allele">"#,
    br#"##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">"#,
    br#"##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period (length of unit)">"#,
    br#"##INFO=<ID=REF,Number=1,Type=Float,Description="Repeat allele length in reference">"#,
];

/// The VCF format lines to be included in the header. See [`make_vcf_header`].
const VCF_FORMAT_LINES: &[&[u8]] = &[
    br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
    br#"##FORMAT=<ID=FT,Number=1,Type=String,Description="Filter tag. Contains PASS if all filters passed, otherwise reason for filter">"#,
    br#"##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">"#,
    br#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of fully spanning reads mapped to locus">"#,
    br#"##FORMAT=<ID=FREQS,Number=1,Type=String,Description="Frequencies observed for each allele length. Keys are allele lengths and values are the number of reads with that allele length.">"#,
    br#"##FORMAT=<ID=REPLEN,Number=1,Type=String,Description="Genotype given in the number of times the unit is repeated for each allele">"#,    
];

/// The VCF filter lines to be included in the header. See [`make_vcf_header`].
const VCF_FILTER_LINES: &[&[u8]] = &[
    br#"##FILTER=<ID=PASS,Description="All filters passed">"#,    
    br#"##FILTER=<ID=UNDEF,Description="Undefined ConSTRain filter">"#,
    br#"##FILTER=<ID=DPZERO,Description="No reads were mapped to locus">"#,
    br#"##FILTER=<ID=DPOOR,Description="Normalised depth of coverage at locus was out of range specified by --min-norm-depth and --max-norm-depth command line arguments">"#,
    br#"##FILTER=<ID=CNZERO,Description="Copy number was zero">"#,
    br#"##FILTER=<ID=CNOOR,Description="Copy number was out of range specified by --max-cn command line argument">"#,
    br#"##FILTER=<ID=CNMISSING,Description="No copy number set for locus. Can happen if contig is missing from karyotype or if only a part of the STR is affected by a CNA">"#,
    br#"##FILTER=<ID=AMBGT,Description="Multiple genotypes are equally likely">"#,
];
