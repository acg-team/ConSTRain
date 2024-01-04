use csv::ReaderBuilder;
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

pub fn trs_to_vcf(tr_regions: &[TandemRepeat]) {

}
