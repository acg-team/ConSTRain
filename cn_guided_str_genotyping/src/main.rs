use ndarray::prelude::*;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::{env, process::exit, collections::HashMap};

use cn_guided_str_genotying::{alelle_lengths_from_cigar, tr_cn_from_cnvs, utils::{self, io_utils}, genotyping::estimate_genotype, repeat};

fn main() {
    // let mut tmp = repeat::TandemRepeat{
    //     reference_info: repeat::RepeatReferenceInfo{seqname: String::from("chr1"), start: 0, end: 1, period: 1, unit: String::from("A")},
    //     copy_number: 5,
    //     allele_lengths: Some(HashMap::from([(9, 1.), (10, 1.), (11, 1.), (12, 1.), (13, 5.), (14, 2.), (15, 23.)])),
    //     genotype: None,
    // };
    // println!("{:?}", tmp);
    // estimate_genotype(&mut tmp);
    // println!("{:?}", tmp);
    // exit(0);
    
    // USAGE: cn-guided-str-genotying <bed file> <ploidy file> <cnv file> <bam file> <n threads>
    // EXAMPLE: cargo run -- ../data/repeats/APC_repeats.bed ../data/genome_architectures/h_sapiens_male.json  ../data/cnv/cnv_test.bed ../data/alignments/duplication_cn3_out1.bam 4
    // EXAMPLE: ./target/debug/cn-guided-str-genotying ../data/repeats/APC_repeats.bed ../data/genome_architectures/h_sapiens_male.json  ../data/cnv/cnv_test.bed ../data/alignments/duplication_cn3_out1.bam 4
    // EXAMPLE: ./target/debug/cn-guided-str-genotying /Users/maxverbiest/PhD/data/str_panels/hg19/tral_and_perf_panel_hg19_reformat.bed ../data/genome_architectures/h_sapiens_male_nochr.json ../data/cnv/cnv_test.bed /Users/maxverbiest/PhD/data/alignments/big_bang_paper/GN.usc_20x.Homo_sapiens_assembly19.fasta.bam 4

    // Eventually make a nicer CLI using Clap crate?
    let args: Vec<String> = env::args().collect();
    let bed_path = &args[1].as_str();
    let ploidy_path = &args[2].as_str();
    let cnv_path = &args[3].as_str();
    let bam_path = &args[4].as_str();
    let n_threads = &args[5].parse::<usize>().unwrap();

    ThreadPoolBuilder::new()
        .num_threads(*n_threads)
        .build_global()
        .unwrap();

    let flanksize = 10;
    let min_reads_per_allele = 10;
    let mut tr_regions = io_utils::trs_from_bed(bed_path, ploidy_path);
    let cnv_regions = io_utils::cnvs_from_bed(cnv_path);

    tr_regions.par_iter_mut().for_each(|tr_region| {
        tr_cn_from_cnvs(tr_region, &cnv_regions);
        if tr_region.copy_number > 0 {
            alelle_lengths_from_cigar(tr_region, bam_path, flanksize);
            estimate_genotype(tr_region, min_reads_per_allele);
        }
        println!("{:?}", tr_region);
    });

    // for r in tr_regions {
    //     println!("{:?}", r);
    // }
}
