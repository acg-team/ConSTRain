use std::{
    collections::HashMap,
    fs, io,
    path::{Path, PathBuf},
    sync::Arc,
};

use constrain::{
    genotyping,
    io::{self as constrain_io, bed::BedFile},
    run,
};
use hex_literal::hex;
use sha2::{Digest, Sha256};

const REPEAT_FILE: &str = "APC_repeats.bed";
const CNV_FILE: &str = "cnv_test.bed";
const BAM_FILE: &str = "duplication_cn3_out1.bam";
const KARYOTYPE_FILE: &str = "h_sapiens_male.json";

fn test_data_dir() -> PathBuf {
    Path::new("tests").join("data")
}

fn sha256_file_digest<P: AsRef<Path>>(path: P) -> Vec<u8> {
    let mut file =
        fs::File::open(&path).expect(&format!("Failed to open file: {}", path.as_ref().display()));
    let mut hasher = Sha256::new();
    _ = io::copy(&mut file, &mut hasher).expect(&format!(
        "Failed to read from file: {}",
        path.as_ref().display()
    ));
    hasher.finalize().to_vec()
}

#[test]
/// Check the input files used for integration tests.
/// If this test fails, it means the repeat file has changed.
/// This is a problem if tests are not updated to reflect the new input file.
fn check_repeat_file() {
    let path = test_data_dir().join(REPEAT_FILE);
    let expect = hex!("3bfac2c7aaebae1103ff48f7bf20318f75bb56f80272e67ac1ebafdada1b5f05");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);
}

#[test]
/// Check the input files used for integration tests.
/// If this test fails, it means the CNV file has changed.
/// This is a problem if tests are not updated to reflect the new input file.
fn check_cnv_file() {
    let path = test_data_dir().join(CNV_FILE);
    let expect = hex!("9020f26e6b5441fa2eb3c4817f0870903672eccbc3cab28fcaa1e9a99e7dbad9");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);
}

#[test]
/// Check the input files used for integration tests.
/// If this test fails, it means the alignment file has changed.
/// This is a problem if tests are not updated to reflect the new input file.
fn check_alignment_file() {
    let path = test_data_dir().join(BAM_FILE);
    let expect = hex!("81c25aace7270c59cc7e56cbf398df14ed8bdf6dcfa33c099f29bb0434909be7");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);
}

#[test]
/// Check the input files used for integration tests.
/// If this test fails, it means the karyotype file has changed.
/// This is a problem if tests are not updated to reflect the new input file.
fn check_karyotype_file() {
    let path = test_data_dir().join(KARYOTYPE_FILE);
    let expect = hex!("61ee3aa08725c85ddba4d63f8b5914486061eb382dbd441948fc96629422898e");
    assert_eq!(sha256_file_digest(path)[..], expect[..]);
}

#[test]
/// Test if the returned allele lengths match the expectations.
fn run_alignment_no_cnv_allele_lengths() {
    let repeat_source = BedFile::new(test_data_dir().join(REPEAT_FILE).to_str().unwrap().into());
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        test_data_dir().join(KARYOTYPE_FILE).to_str().unwrap(),
        20,
        None::<&BedFile>,
    )
    .unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(
        &mut tr_regions,
        &partitions_map,
        test_data_dir().join(BAM_FILE).to_str().unwrap(),
        None,
        5,
        1.0,
        None,
        0,
    )
    .unwrap();

    for repeat in tr_regions.iter() {
        println!("{:?}", repeat);
    }

    let expect_allele_lengths = vec![
        Some(HashMap::from([(10, 44.0), (15, 39.0), (16, 41.0)])),
        Some(HashMap::from([(9, 145.0)])),
        Some(HashMap::from([(18, 95.0)])),
        Some(HashMap::from([(26, 109.0)])),
        Some(HashMap::from([(14, 145.0)])),
    ];

    for (repeat, expect) in tr_regions.iter().zip(expect_allele_lengths) {
        assert_eq!(repeat.allele_lengths, expect);
    }
}

#[test]
/// Test if the returned genotypes match the expectations.
fn run_alignment_no_cnv_genotypes() {
    let repeat_source = BedFile::new(test_data_dir().join(REPEAT_FILE).to_str().unwrap().into());
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        test_data_dir().join(KARYOTYPE_FILE).to_str().unwrap(),
        20,
        None::<&BedFile>,
    )
    .unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(
        &mut tr_regions,
        &partitions_map,
        test_data_dir().join(BAM_FILE).to_str().unwrap(),
        None,
        5,
        1.0,
        None,
        0,
    )
    .unwrap();

    let expect_genotypes = vec![
        Some(vec![(10, 1.0), (16, 1.0)]),
        Some(vec![(9, 2.0)]),
        Some(vec![(18, 2.0)]),
        Some(vec![(26, 2.0)]),
        Some(vec![(14, 2.0)]),
    ];

    for (repeat, expect) in tr_regions.iter().zip(expect_genotypes) {
        assert_eq!(repeat.genotype, expect);
    }
}

#[test]
/// Test if the returned filter tags match the expectations.
fn run_alignment_no_cnv_filter_tags() {
    let repeat_source = BedFile::new(test_data_dir().join(REPEAT_FILE).to_str().unwrap().into());
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        test_data_dir().join(KARYOTYPE_FILE).to_str().unwrap(),
        20,
        None::<&BedFile>,
    )
    .unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(
        &mut tr_regions,
        &partitions_map,
        test_data_dir().join(BAM_FILE).to_str().unwrap(),
        None,
        5,
        1.0,
        None,
        0,
    )
    .unwrap();

    let expect_tags = vec!["PASS", "PASS", "PASS", "PASS", "PASS"];

    for (repeat, expect) in tr_regions.into_iter().zip(expect_tags) {
        assert_eq!(repeat.filter.name(), expect);
    }
}

#[test]
/// Test if the returned allele lengths match the expectations.
fn run_alignment_cnv_allele_lengths() {
    let repeat_source = BedFile::new(test_data_dir().join(REPEAT_FILE).to_str().unwrap().into());
    let cnv_source = BedFile::new(test_data_dir().join(CNV_FILE).to_str().unwrap().into());
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        test_data_dir().join(KARYOTYPE_FILE).to_str().unwrap(),
        20,
        Some(&cnv_source),
    )
    .unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(
        &mut tr_regions,
        &partitions_map,
        test_data_dir().join(BAM_FILE).to_str().unwrap(),
        None,
        5,
        1.0,
        None,
        0,
    )
    .unwrap();

    let expect_allele_lengths = vec![
        Some(HashMap::from([(10, 44.0), (15, 39.0), (16, 41.0)])),
        None,
        Some(HashMap::from([(18, 95.0)])),
        None,
        Some(HashMap::from([(14, 145.0)])),
    ];

    for (repeat, expect) in tr_regions.iter().zip(expect_allele_lengths) {
        assert_eq!(repeat.allele_lengths, expect);
    }
}

#[test]
/// Test if the returned genotypes match the expectations.
fn run_alignment_cnv_genotypes() {
    let repeat_source = BedFile::new(test_data_dir().join(REPEAT_FILE).to_str().unwrap().into());
    let cnv_source = BedFile::new(test_data_dir().join(CNV_FILE).to_str().unwrap().into());
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        test_data_dir().join(KARYOTYPE_FILE).to_str().unwrap(),
        20,
        Some(&cnv_source),
    )
    .unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(
        &mut tr_regions,
        &partitions_map,
        test_data_dir().join(BAM_FILE).to_str().unwrap(),
        None,
        5,
        1.0,
        None,
        0,
    )
    .unwrap();

    let expect_genotypes = vec![
        Some(vec![(10, 1.0), (15, 1.0), (16, 1.0)]),
        None,
        Some(vec![(18, 1.0)]),
        None,
        Some(vec![(14, 2.0)]),
    ];

    for (repeat, expect) in tr_regions.iter().zip(expect_genotypes) {
        assert_eq!(repeat.genotype, expect);
    }
}

#[test]
/// Test if the returned filter tags match the expectations.
fn run_alignment_cnv_filter_tags() {
    let repeat_source = BedFile::new(test_data_dir().join(REPEAT_FILE).to_str().unwrap().into());
    let cnv_source = BedFile::new(test_data_dir().join(CNV_FILE).to_str().unwrap().into());
    let (mut tr_regions, observed_copy_numbers) = constrain_io::load_tandem_repeats(
        &repeat_source,
        test_data_dir().join(KARYOTYPE_FILE).to_str().unwrap(),
        20,
        Some(&cnv_source),
    )
    .unwrap();

    let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

    run(
        &mut tr_regions,
        &partitions_map,
        test_data_dir().join(BAM_FILE).to_str().unwrap(),
        None,
        5,
        1.0,
        None,
        0,
    )
    .unwrap();

    let expect_tags = vec!["PASS", "CNMISSING", "PASS", "CNMISSING", "PASS"];

    for (repeat, expect) in tr_regions.into_iter().zip(expect_tags) {
        assert_eq!(repeat.filter.name(), expect);
    }
}
