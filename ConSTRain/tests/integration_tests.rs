use std::{fs, io, path::Path, sync::Arc};

use hex_literal::hex;
use sha2::{Sha256, Digest};
// cargo run --bin ConSTRain -- alignment --repeats ../data/repeats/APC_repeats.bed --karyotype ../data/genome_architectures/h_sapiens_male.json --alignment ../data/alignments/duplication_cn3_out1.bam --cnvs ../data/cnv/cnv_test.bed
use constrain::{genotyping, io::{self as constrain_io, bed::BedFile}, run};

const TEST_DATA_DIR: &str = "./tests/data/";
const REPEAT_FILE: &str = "APC_repeats.bed";
const CNV_FILE: &str = "cnv_test.bed";
const BAM_FILE: &str = "duplication_cn3_out1.bam";
const KARYOTYPE_FILE: &str = "h_sapiens_male.json";

fn sha256_file_digest<P: AsRef<Path>>(path: P) -> Vec<u8> {
    let mut file = fs::File::open(&path).expect(&format!("Failed to open file: {}", path.as_ref().display()));
    let mut hasher = Sha256::new();
    _ = io::copy(&mut file, &mut hasher).expect(&format!("Failed to read from file: {}", path.as_ref().display()));
    hasher.finalize().to_vec()
}

#[test]
/// Check the input files used for integration tests.
/// If this test fails, it means one or more of the input files have changed.
/// This is a problem if tests are not updated to reflect the new input files.
fn check_input_files() {
    // Check bed file with repeats
    let path = Path::new(TEST_DATA_DIR).join(REPEAT_FILE);
    let expect = hex!("bf209733f318ed6b79b62dce2287f4ee3049fae47a36c7e56c0c72c6bb3e6db9");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);

    // Check bed file with CNVs
    let path = Path::new(TEST_DATA_DIR).join(CNV_FILE);
    let expect = hex!("9020f26e6b5441fa2eb3c4817f0870903672eccbc3cab28fcaa1e9a99e7dbad9");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);

    // Check BAM alignment file
    let path = Path::new(TEST_DATA_DIR).join(BAM_FILE);
    let expect = hex!("81c25aace7270c59cc7e56cbf398df14ed8bdf6dcfa33c099f29bb0434909be7");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);

    // Check karyotype file
    let path = Path::new(TEST_DATA_DIR).join(KARYOTYPE_FILE);
    let expect = hex!("61ee3aa08725c85ddba4d63f8b5914486061eb382dbd441948fc96629422898e");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);
}

#[test]
fn run_alignment_no_cnv() {
    let repeat_source = BedFile::new(format!("{TEST_DATA_DIR}/{REPEAT_FILE}"));
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        &format!("{TEST_DATA_DIR}/{KARYOTYPE_FILE}"),
        20,
        None::<&BedFile>,
    ).unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(&mut tr_regions, &partitions_map, &format!("{TEST_DATA_DIR}/{BAM_FILE}"), None, 5, 1.0, None, 0).unwrap();
}

