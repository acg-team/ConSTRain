use rayon::{prelude::*, ThreadPoolBuilder, current_thread_index};
use rust_htslib::{bam::IndexedReader};
use std::{env, sync::Arc, thread, time};

use cn_guided_str_genotying::{
    fetch_allele_lengths, genotyping::estimate_genotype, make_compositions_map, tr_cn_from_cnvs,
    utils::io_utils,
};

fn main() {
    // USAGE: cn-guided-str-genotying <bed file> <ploidy file> <cnv file> <bam file> <n threads>
    // EXAMPLE: cargo run -- ../data/repeats/APC_repeats.bed ../data/genome_architectures/h_sapiens_male.json  ../data/cnv/cnv_test.bed ../data/alignments/duplication_cn3_out1.bam 4
    // EXAMPLE: ./target/debug/cn-guided-str-genotying ../data/repeats/APC_repeats.bed ../data/genome_architectures/h_sapiens_male.json  ../data/cnv/cnv_test.bed ../data/alignments/duplication_cn3_out1.bam 4
    // EXAMPLE: ./target/release/cn-guided-str-genotying /Users/maxverbiest/PhD/data/str_panels/hg19/tral_and_perf_panel_hg19_reformat.bed ../data/genome_architectures/h_sapiens_male_nochr.json ../data/cnv/cnv_test.bed /Users/maxverbiest/PhD/data/alignments/big_bang_paper/GN.usc_20x.Homo_sapiens_assembly19.fasta.bam 4
    // EXAMPLE: ./target/release/cn-guided-str-genotying /Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_updated.bed ../data/genome_architectures/h_sapiens_female.json ../data/cnv/cnv_test.bed  /Users/maxverbiest/PhD/data/alignments/1000g/HG00138.final.cram 4

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

    // Currently, io_utils functions return all copy numbers they encounter
    // while reading data. This is then used to create compositions_map.
    // Could instead switch to Arc<RwLock<HashMap<..., ...>>> implementation
    // where compositions_map is updated if a new copy number is encountered.
    let (copy_numbers, mut tr_regions) = io_utils::trs_from_bed(bed_path, ploidy_path);
    eprintln!("Read {} TR regions", tr_regions.len());
    let (cnv_copy_numbers, cnv_regions) = io_utils::cnvs_from_bed(cnv_path);
    eprintln!("Read {} CNVs", cnv_regions.len());
    let copy_numbers: Vec<usize> = copy_numbers.union(&cnv_copy_numbers).filter_map(|x| if *x < 20 {Some(*x)} else {None}).collect();
    eprintln!("Generating compositions for copy numbers {:?}", copy_numbers);
    let compositions_map = Arc::new(make_compositions_map(&copy_numbers));
    eprintln!("Generated compositions");

    let reference = "/Users/maxverbiest/PhD/data/genomes/GRCh38.d1.vd1.fa";    

    let chunksize = tr_regions.len() / n_threads + 1;
    tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
        let mut bam = IndexedReader::from_path(bam_path).unwrap();
        bam.set_reference(reference).unwrap();
        
        for tr_region in tr_regions {
            tr_cn_from_cnvs(tr_region, &cnv_regions);
            let fetch_request = tr_region.reference_info.get_fetch_definition();     

            let fetch_result = bam.fetch(fetch_request);
            if fetch_result.is_err() {
                eprintln!("Fetching returned error for TR: {:?}", tr_region);
                continue;
            }

            fetch_allele_lengths(tr_region, &mut bam, flanksize);
            estimate_genotype(
                tr_region,
                min_reads_per_allele,
                Arc::clone(&compositions_map),
            );
            println!("{:?}", tr_region);
        }
        eprintln!("Finished chunk on thread {}", current_thread_index().unwrap());
        thread::sleep(time::Duration::from_secs(1));
        // Segmentation fault after this point when reading from CRAM file
        // I guess it's because of something that happens when the IndexedReader goes
        // out of scope and gets dropped
        eprintln!("Attempting to drop IndexedReader on thread {}", current_thread_index().unwrap());        
        drop(bam);
        eprintln!("Dropped IndexedReader on thread {}", current_thread_index().unwrap());        
    });

    // for r in tr_regions {
    //     println!("{:?}", r);
    // }
}
