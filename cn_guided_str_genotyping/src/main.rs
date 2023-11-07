use rayon::ThreadPoolBuilder;
use std::env;

use cn_guided_str_genotying;


fn main() {
    // USAGE: cn-guided-str-genotying <bed file> <bam file> <threads>
    // EXAMPLE: ./target/debug/cn-guided-str-genotying ../data/repeats/APC_repeats.bed ../data/alignments/patient_3.bam ../data/genome_architectures/h_sapiens_male.json 4
    // EXAMPLE: cargo run -- ../data/repeats/APC_repeats.tsv ../data/alignments/patient_3.bam ../data/genome_architectures/h_sapiens_male.json 4
    // Eventually make a nicer CLI using Clap crate?
    let args: Vec<String> = env::args().collect();
    let bed_path = &args[1];
    let bam_path = &args[2];
    let ploidy_path = &args[3];
    let n_threads = &args[4].parse::<usize>().unwrap();

    ThreadPoolBuilder::new().num_threads(*n_threads).build_global().unwrap();

    let flanksize = 4;
    cn_guided_str_genotying::call_tr_allele_lengths(bed_path, bam_path, ploidy_path, flanksize);
}
