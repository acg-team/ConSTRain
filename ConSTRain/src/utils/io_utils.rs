//! # IO Utils
//!
//! Home of I/O functionality needed in the ConSTRain library. Provides
//! functionality to serialize/deserialze tandem repeats from BED files,
//! VCF files.
use csv::ReaderBuilder;
use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Record, Writer};
use serde_json::Value;
use std::io;
use std::str::Utf8Error;
use std::{collections::HashSet, error::Error, fs::File, io::BufReader};

use crate::repeat::{RepeatReferenceInfo, TandemRepeat};
use crate::utils::CopyNumberVariant;

/// Read tandem repeat regions specified in the bed file at `bed_path` into `tr_buffer`.
/// Copy numbers of tandem repeats are set based on the values specified in `ploidy_path`, which should be a
/// json file mapping contig identifiers to copy number values. All copy number values that are observed while
/// reading TR regions are added to `cn_buffer`, which can be useful later to generate partitions for only the
/// relevant (i.e., observed) copy numbers.
pub fn trs_from_bed(
    bed_path: &str,
    ploidy_path: &str,
    tr_buffer: &mut Vec<TandemRepeat>,
    cn_buffer: &mut HashSet<usize>,
) -> Result<(), Box<dyn Error>> {
    let ploidy = ploidy_from_json(ploidy_path)?;

    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(bed_path)?;

    for result in bed_reader.deserialize() {
        let ref_info: RepeatReferenceInfo = result?;

        // Get base copy number of contig from ploidy json.
        // If the contig does not exist in the json, default to CN 0 (which means TRs on this contig won't be genotyped)
        let cn = ploidy[&ref_info.seqname].as_u64().unwrap_or(0) as usize;
        if cn == 0 {
            eprintln!("WARNING: contig for TR '{}' was not found in the ploidy json file. Setting CN to 0", ref_info.get_fetch_definition_s())
        }

        let tr_region = TandemRepeat {
            reference_info: ref_info,
            copy_number: cn,
            allele_lengths: None,
            genotype: None,
        };
        cn_buffer.insert(cn);
        tr_buffer.push(tr_region);
    }

    Ok(())
}

/// Read contig-level baseline copy number values from a json file at `path`.
/// The json should contain contig names as keys and integer copy numbers as values, e.g.:
/// ```
/// {
///     "chr1": 2,
///     ... other chromosomes ...
///     "chrY": 0
/// }
/// ```
fn ploidy_from_json(path: &str) -> Result<Value, io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    Ok(serde_json::from_reader(reader)?)
}

/// Read copy number variants specified in the bed file at `cnv_path` into `cnv_buffer`.
/// All copy number values that are observed while reading TR regions are added to `cn_buffer`,
/// which can be useful later to generate partitions for only the relevant (i.e., observed) copy numbers.
pub fn cnvs_from_bed(
    cnv_path: &str,
    cnv_buffer: &mut Vec<CopyNumberVariant>,
    cn_buffer: &mut HashSet<usize>,
) -> Result<(), Box<dyn Error>> {
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(cnv_path)?;

    for result in bed_reader.deserialize() {
        let cnv: CopyNumberVariant = result?;
        cn_buffer.insert(cnv.cn);
        cnv_buffer.push(cnv);
    }

    Ok(())
}

/// Write (genotyped) tandem repeat regions to a vcf file.
pub fn trs_to_vcf(
    tr_regions: &[TandemRepeat],
    targets: &[&[u8]],
    lengths: &[u64],
    sample_name: &str,
) -> Result<(), Box<dyn Error>> {
    // Create header, write uncompressed VCF to stdout (header gets written when Writer is constructed)
    let header = make_vcf_header(targets, lengths, sample_name)?;
    let mut vcf = Writer::from_stdout(&header, true, Format::Vcf)?;

    for tr_region in tr_regions {
        // Create record for repeat and add reference information
        let mut record = vcf.empty_record();
        add_reference_info(&mut record, &vcf, tr_region)?;

        // Using the estimated genotype for this locus, create the alleles and genotype string to add to the VCF
        add_alleles_genotypes(&mut record, tr_region)?;

        // Add additional information that was extracted from the alignment to the VCF af FORMAT fields
        add_additional_info(&mut record, tr_region)?;

        // Write record
        vcf.write(&record)?
    }

    Ok(())
}

/// The VCF info lines to be included in the header. See [make_vcf_header].
const VCF_INFO_LINES: &[&[u8]] = &[
    br#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of reference allele">"#,
    br#"##INFO=<ID=RU,Number=1,Type=String,Description="Repeat motif">"#,
    br#"##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period (length of motif)">"#,
    br#"##INFO=<ID=REF,Number=1,Type=Float,Description="Repeat allele length in reference">"#,
];

/// The VCF format lines to be included in the header. See [make_vcf_header].
const VCF_FORMAT_LINES: &[&[u8]] = &[
    br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
    br#"##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">"#,
    br#"##FORMAT=<ID=FREQS,Number=1,Type=String,Description="Frequencies observed for each allele length. Keys are allele lengths and values are the number of reads with that allele length.">"#,
    br#"##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Genotype given in number of copies of the repeat motif">"#,
];

/// Construct VCF a header. First, include information about the target contigs, followed by the [VCF_INFO_LINES] and [VCF_FORMAT_LINES].
/// Then, write TR variant calls to stdout.
fn make_vcf_header(
    targets: &[&[u8]],
    lengths: &[u64],
    sample_name: &str,
) -> Result<Header, Utf8Error> {
    let mut header = Header::new();

    for (target, length) in targets.iter().zip(lengths.iter()) {
        // let target_str = std::str::from_utf8(*target).expect("Error parsing contig name!");
        let target_str = std::str::from_utf8(target)?;
        let header_contig_line = format!(r#"##contig=<ID={},length={}>"#, target_str, length);
        header.push_record(header_contig_line.as_bytes());
    }

    for header_info_line in VCF_INFO_LINES.iter() {
        header.push_record(header_info_line);
    }
    for header_format_line in VCF_FORMAT_LINES.iter() {
        header.push_record(header_format_line);
    }
    header.push_sample(sample_name.as_bytes());

    Ok(header)
}

/// Add information about how the `tr_region` looks in the reference genome to the
/// VCF `record`. We need the `vcf` struct to get the contig id from the VCF header.
fn add_reference_info(
    record: &mut Record,
    vcf: &Writer,
    tr_region: &TandemRepeat,
) -> Result<(), Box<dyn Error>> {
    let rid = vcf
        .header()
        .name2rid(tr_region.reference_info.seqname.as_bytes())?;
    let ref_len = tr_region.reference_info.get_reference_len() as usize;
    record.set_rid(Some(rid));
    record.set_pos(tr_region.reference_info.start);

    record.push_info_integer(b"END", &[tr_region.reference_info.end as i32])?;
    record.push_info_string(b"RU", &[tr_region.reference_info.unit.as_bytes()])?;
    record.push_info_integer(b"PERIOD", &[tr_region.reference_info.period as i32])?;
    record.push_info_integer(b"REF", &[ref_len as i32])?;

    Ok(())
}

/// Add genotyping results stored on `tr_region` to the VCF `record`.
fn add_alleles_genotypes(
    record: &mut Record,
    tr_region: &TandemRepeat,
) -> Result<(), Box<dyn Error>> {
    let ref_len = tr_region.reference_info.get_reference_len() as usize;
    let ref_allele = tr_region.reference_info.unit.repeat(ref_len);
    let mut alleles = vec![ref_allele];
    let mut genotype: Vec<GenotypeAllele> = Vec::new();
    match &tr_region.genotype {
        None => genotype.push(GenotypeAllele::UnphasedMissing), // No genotype estimated for repeat
        Some(gts) => {
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
    let alleles = Vec::from_iter(alleles.iter().map(|i| i.as_bytes()));
    record.set_alleles(&alleles)?;
    record.push_genotypes(&genotype)?;
    Ok(())
}

/// Add additional information stored on `tr_region` to the VCF `record`
fn add_additional_info(
    record: &mut Record,
    tr_region: &TandemRepeat,
) -> Result<(), Box<dyn Error>> {
    record.push_format_integer(b"CN", &[tr_region.copy_number as i32])?;

    let mut allele_freqs = tr_region.allele_freqs_as_tuples();
    allele_freqs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let allele_freqs: Vec<String> = allele_freqs
        .iter()
        .map(|i| format!("{},{}", i.0, i.1))
        .collect();
    let allele_freqs = allele_freqs.join("|");
    record.push_format_string(b"FREQS", &[allele_freqs.as_bytes()])?;

    let gt_as_allele_lens: Vec<String> = tr_region
        .gt_as_allele_lengths()
        .iter()
        .map(|i| i.to_string())
        .collect();
    let gt_as_allele_lens = gt_as_allele_lens.join(",");
    record.push_format_string(b"REPCN", &[gt_as_allele_lens.as_bytes()])?;

    Ok(())
}

#[cfg(test)]
mod tests {
    // use super::*;

    // #[test]
    // fn adding_alleles_genotypes() {
    //     let ref_info = RepeatReferenceInfo{
    //         seqname: String::from("chr_test"),
    //         start: 10,
    //         end: 20,
    //         period: 2,
    //         unit: String::from("AT")
    //     };
    //     let tr_region = TandemRepeat{
    //         reference_info: ref_info,
    //         copy_number: 2,
    //         allele_lengths: Some(HashMap::from([
    //             (9, 20.),
    //             (10, 20.),
    //             (11, 20.),
    //         ])),
    //         genotype: Some(vec![(9, 1.), (10, 1.), (11, 1.)])
    //     };
    // }
}
