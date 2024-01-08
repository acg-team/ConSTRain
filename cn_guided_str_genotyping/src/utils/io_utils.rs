use csv::ReaderBuilder;
use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Record, Writer};
use serde_json::Value;
use std::{collections::HashSet, fs::File, io::BufReader};

use crate::repeat::{RepeatReferenceInfo, TandemRepeat};
use crate::utils::CopyNumberVariant;

pub fn trs_from_bed(bed_path: &str, ploidy_path: &str) -> (HashSet<usize>, Vec<TandemRepeat>) {
    let ploidy = ploidy_from_json(ploidy_path);

    let mut tr_regions = Vec::new();
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(bed_path)
        .unwrap();

    let mut copy_numbers: HashSet<usize> = HashSet::new();
    for result in bed_reader.deserialize() {
        let ref_info: RepeatReferenceInfo = result.unwrap();
        let cn = ploidy[ref_info.seqname.as_str()].as_u64().unwrap_or(0) as usize;

        let tr_region = TandemRepeat {
            reference_info: ref_info,
            copy_number: cn,
            allele_lengths: None,
            genotype: None,
        };

        copy_numbers.insert(cn);
        tr_regions.push(tr_region);
    }

    (copy_numbers, tr_regions)
}

fn ploidy_from_json(path: &str) -> Value {
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let v: Value = serde_json::from_reader(reader).unwrap();

    v
}

pub fn cnvs_from_bed(cnv_path: &str) -> (HashSet<usize>, Vec<CopyNumberVariant>) {
    let mut cnv_regions = Vec::new();
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(cnv_path)
        .unwrap();

    let mut copy_numbers: HashSet<usize> = HashSet::new();
    for result in bed_reader.deserialize() {
        let cnv: CopyNumberVariant = result.unwrap();
        copy_numbers.insert(cnv.cn);
        cnv_regions.push(cnv);
    }

    (copy_numbers, cnv_regions)
}

pub fn trs_to_vcf(
    tr_regions: &[TandemRepeat],
    targets: &[&[u8]],
    lengths: &[u64],
    sample_name: &str,
) {
    // Create header, write uncompressed VCF to stdout (header gets written when Writer is constructed)
    let header = make_vcf_header(targets, lengths, sample_name);
    let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();

    for tr_region in tr_regions {
        // Create record for repeat and add reference information
        let mut record = vcf.empty_record();
        add_reference_info(&mut record, &vcf, tr_region);

        // Using the estimated genotype for this locus, create the alleles and genotype string to add to the VCF
        add_alleles_genotypes(&mut record, tr_region);

        // Add additional information that was extracted from the alignment to the VCF af FORMAT fields
        add_additional_info(&mut record, tr_region);

        // Write record
        vcf.write(&record).unwrap()
    }
}

const VCF_INFO_LINES: &'static [&[u8]] = &[
    br#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of reference allele">"#,
    br#"##INFO=<ID=RU,Number=1,Type=String,Description="Repeat motif">"#,
    br#"##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period (length of motif)">"#,
    br#"##INFO=<ID=REF,Number=1,Type=Float,Description="Repeat allele length in reference">"#,
];

const VCF_FORMAT_LINES: &'static [&[u8]] = &[
    br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
    br#"##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">"#,
    br#"##FORMAT=<ID=FREQS,Number=1,Type=String,Description="Frequencies observed for each allele length. Keys are allele lengths and values are the number of reads with that allele length.">"#,
];

fn make_vcf_header(targets: &[&[u8]], lengths: &[u64], sample_name: &str) -> Header {
    // Create VCF header, populate with contigs, INFO, FORMAT lines, and sample name
    let mut header = Header::new();

    for (target, length) in targets.iter().zip(lengths.iter()) {
        let target_str = std::str::from_utf8(*target).expect("Error parsing contig name!");
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

    header
}

fn add_reference_info(record: &mut Record, vcf: &Writer, tr_region: &TandemRepeat) {
    let rid = vcf
        .header()
        .name2rid(tr_region.reference_info.seqname.as_bytes())
        .expect("Failed to set reference ID");
    let ref_len = tr_region.reference_info.get_reference_len() as usize;
    record.set_rid(Some(rid));
    record.set_pos(tr_region.reference_info.start as i64);

    record
        .push_info_integer(b"END", &[tr_region.reference_info.end as i32])
        .expect("Failed to set END format field");
    record
        .push_info_string(b"RU", &[tr_region.reference_info.unit.as_bytes()])
        .expect("Failed to set RU format field");
    record
        .push_info_integer(b"PERIOD", &[tr_region.reference_info.period as i32])
        .expect("Failed to set PERIOD format field");
    record
        .push_info_integer(b"REF", &[ref_len as i32])
        .expect("Failed to set REF format field");
}

fn add_alleles_genotypes(record: &mut Record, tr_region: &TandemRepeat) {
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
    record.set_alleles(&alleles).expect("Failed to set alleles");
    record
        .push_genotypes(&genotype)
        .expect("Failed to set genotype");
}

fn add_additional_info(record: &mut Record, tr_region: &TandemRepeat) {
    record
        .push_format_integer(b"CN", &[tr_region.copy_number as i32])
        .expect("Failed to set copy number");
    let allele_freqs: Vec<String> = tr_region
        .allele_freqs_as_tuples()
        .iter()
        .map(|i| format!("{},{}", i.0, i.1))
        .collect();
    let allele_freqs = allele_freqs.join("|");
    record
        .push_format_string(b"FREQS", &[allele_freqs.as_bytes()])
        .expect("Failed to set allele frequencies");
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
