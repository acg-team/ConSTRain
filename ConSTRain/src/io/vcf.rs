use std::{
    collections::{HashMap, HashSet},
    str,
};

use anyhow::{bail, Context, Result};
use log::info;
use rust_htslib::bcf::{
    header::{Header, HeaderView},
    record::GenotypeAllele,
    Format, Read, Reader, Record, Writer,
};

use crate::{
    karyotype::Karyotype,
    repeat::{RepeatReferenceInfo, TandemRepeat},
    utils::VcfFilter,
};

/// Read tandem repeat regions specified in the vcf file at `vcf_path` into `tr_buffer`.
/// Copy numbers of tandem repeats are set based on the values specified in `karyotype`. All copy number
/// values that are observed while reading TR regions are added to `observed_cn_buffer`, which is used later to
/// generate partitions for only the relevant (i.e., observed) copy numbers.
pub fn read_trs(
    vcf_path: &str,
    karyotype: &Karyotype,
    tr_buffer: &mut Vec<TandemRepeat>,
    observed_cn_buffer: &mut HashSet<usize>,
    sample_name: &str,
) -> Result<()> {
    let mut bcf = Reader::from_path(vcf_path)
        .with_context(|| format!("Failed to open VCF file at {vcf_path}"))?;

    let header = bcf.header().to_owned();
    let sample_idx = header
        .sample_id(sample_name.as_bytes())
        .with_context(|| format!("Sample {sample_name} not found in file {vcf_path}"))?;

    for record in bcf.records() {
        let record =
            record.with_context(|| format!("Error reading VCF record in file {vcf_path}"))?;
        let ref_info = ref_info_from_record(&record, &header)?;

        // Get allele lengths from Record, parse into hashmap
        let allele_freqs = allele_lens_from_record(&record, sample_idx)?;

        // Get GT from Record, parse into vector
        // let genotype = genotype_from_record(&record, sample_idx)?;
        let genotype = None;

        let mut tr = TandemRepeat {
            reference_info: ref_info,
            copy_number: 0, // placeholder copy number value
            allele_lengths: allele_freqs,
            genotype: genotype,
            filter: VcfFilter::Pass,
        };
        tr.set_cn_from_karyotpe(karyotype);
        observed_cn_buffer.insert(tr.copy_number);
        tr_buffer.push(tr);
    }

    info!("Read {} TR regions from {vcf_path}", tr_buffer.len());
    Ok(())
}

/// Is this really the way to get info fields from a record??? Seems incredibly tedious
fn get_info_int(record: &Record, tag: &str) -> Result<i64> {
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

fn get_info_str(record: &Record, tag: &str) -> Result<String> {
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

/// Might be defined on the RepeatReferenceInfo struct instead
fn ref_info_from_record(record: &Record, header: &HeaderView) -> Result<RepeatReferenceInfo> {
    let rid = record.rid().context("Failed to get record rid")?;
    let contig = str::from_utf8(header.rid2name(rid)?)?.to_string();
    let end = get_info_int(record, "END")?;
    let period = get_info_int(record, "PERIOD")?;
    let unit = get_info_str(record, "RU")?;

    Ok(RepeatReferenceInfo {
        seqname: contig,
        start: record.pos(),
        end,
        period,
        unit,
    })
}

fn allele_lens_from_record(
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

/// Write genotyped tandem repeat regions to a vcf file, reusing header from vcf file at `vcf_path`.
pub fn write_reuse_header(tr_regions: &[TandemRepeat], vcf_path: &str) -> Result<()> {
    let reader = Reader::from_path(vcf_path)
        .with_context(|| format!("Failed to open VCF file at {vcf_path}"))?;
    let header = Header::from_template(reader.header());
    let mut vcf = Writer::from_stdout(&header, true, Format::Vcf)?;

    for tr_region in tr_regions {
        // Create record for repeat and add reference information
        let mut record = vcf.empty_record();
        add_info_fields(&mut record, &vcf, tr_region)?;

        // Using the estimated genotype for this locus, create the alleles and genotype string to add to the VCF
        add_alleles_genotypes(&mut record, tr_region)?;

        // Add additional information that was extracted from the alignment to the VCF af FORMAT fields
        add_format_fields(&mut record, tr_region)?;

        // Write record
        vcf.write(&record)?;
    }

    Ok(())
}

/// Write genotyped tandem repeat regions to a vcf file.
pub fn write(
    tr_regions: &[TandemRepeat],
    targets: &[String],
    lengths: &[u64],
    sample_name: &str,
) -> Result<()> {
    // Create header, write uncompressed VCF to stdout (header gets written when Writer is constructed)
    let header = make_vcf_header(targets, lengths, sample_name);
    let mut vcf = Writer::from_stdout(&header, true, Format::Vcf)?;

    for tr_region in tr_regions {
        // Create record for repeat and add reference information
        let mut record = vcf.empty_record();
        add_info_fields(&mut record, &vcf, tr_region)?;

        // Using the estimated genotype for this locus, create the alleles and genotype string to add to the VCF
        add_alleles_genotypes(&mut record, tr_region)?;

        // Add additional information that was extracted from the alignment to the VCF af FORMAT fields
        add_format_fields(&mut record, tr_region)?;

        // Write record
        vcf.write(&record)?;
    }

    Ok(())
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

/// Construct VCF a header. First, include information about the target contigs, followed by the [`VCF_INFO_LINES`], [`VCF_FORMAT_LINES`].
/// Then, write TR variant calls to stdout.
fn make_vcf_header(targets: &[String], lengths: &[u64], sample_name: &str) -> Header {
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

/// Add information about how the `tr_region` looks in the reference genome to the
/// VCF `record`. We need the `vcf` struct to get the contig id from the VCF header.
fn add_info_fields(record: &mut Record, vcf: &Writer, tr_region: &TandemRepeat) -> Result<()> {
    let context = || {
        format!(
            "Error setting INFO field value for {}",
            tr_region.reference_info.get_fetch_definition_s()
        )
    };

    let rid = vcf
        .header()
        .name2rid(tr_region.reference_info.seqname.as_bytes())
        .with_context(context)?;
    let ref_len = tr_region.reference_info.get_reference_len() as usize;
    record.set_rid(Some(rid));
    record.set_pos(tr_region.reference_info.start);

    record
        .push_info_integer(b"END", &[tr_region.reference_info.end as i32])
        .with_context(context)?;
    record
        .push_info_string(b"RU", &[tr_region.reference_info.unit.as_bytes()])
        .with_context(context)?;
    record
        .push_info_integer(b"PERIOD", &[tr_region.reference_info.period as i32])
        .with_context(context)?;
    record
        .push_info_integer(b"REF", &[ref_len as i32])
        .with_context(context)?;

    Ok(())
}

/// Add genotyping results stored on `tr_region` to the VCF `record`.
fn add_alleles_genotypes(record: &mut Record, tr_region: &TandemRepeat) -> Result<()> {
    let ref_len = tr_region.reference_info.get_reference_len() as usize;
    let ref_allele = tr_region.reference_info.unit.repeat(ref_len);
    let mut alleles = vec![ref_allele];
    let mut genotype: Vec<GenotypeAllele> = Vec::new();

    match &tr_region.genotype {
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
                let allele = tr_region.reference_info.unit.repeat(len as usize);
                alleles.push(allele);
            }
        }
    }

    let alleles: Vec<&[u8]> = alleles.iter().map(std::string::String::as_bytes).collect();
    record.set_alleles(&alleles).with_context(|| {
        format!(
            "Error setting alleles for {}",
            tr_region.reference_info.get_fetch_definition_s()
        )
    })?;
    record.push_genotypes(&genotype).with_context(|| {
        format!(
            "Error setting genotype for {}",
            tr_region.reference_info.get_fetch_definition_s()
        )
    })?;
    Ok(())
}

/// Add additional information stored on `tr_region` to the VCF `record` format fields
fn add_format_fields(record: &mut Record, tr_region: &TandemRepeat) -> Result<()> {
    let context = || {
        format!(
            "Error setting FORMAT field value for {}",
            tr_region.reference_info.get_fetch_definition_s()
        )
    };

    record
        .push_format_string(b"FT", &[tr_region.filter.name().as_bytes()])
        .with_context(context)?;

    if !matches!(tr_region.filter, VcfFilter::CnMissing) {
        record
            .push_format_integer(b"CN", &[tr_region.copy_number as i32])
            .with_context(context)?;
    }
    
    if let Some(depth) = tr_region.get_n_mapped_reads() {
        record.push_format_integer(b"DP", &[depth as i32]).with_context(context)?;
    } else {
        record.push_format_integer(b"DP", &[0]).with_context(context)?;
    }

    let mut allele_freqs = tr_region.allele_freqs_as_tuples();
    allele_freqs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let allele_freqs: Vec<String> = allele_freqs
        .iter()
        .map(|i| format!("{},{}", i.0, i.1))
        .collect();
    let allele_freqs = allele_freqs.join("|");
    record
        .push_format_string(b"FREQS", &[allele_freqs.as_bytes()])
        .with_context(context)?;

    let gt_as_allele_lens: Vec<String> = tr_region
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
