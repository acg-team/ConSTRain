use rayon::{current_thread_index, prelude::*, ThreadPoolBuilder};
use rust_htslib::{bam, htslib, utils};
use std::path::Path;
use std::{env, ffi, sync::Arc};

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
    // EXAMPLE: ./target/release/cn-guided-str-genotying ../data/repeats/HG00138_strs_in_cnv.bed ../data/genome_architectures/h_sapiens_male.json ../data/cnv/HG00138_pgx.bed  /Users/maxverbiest/PhD/data/alignments/1000g/HG00138.final.cram 8

    // Eventually make a nicer CLI using Clap crate?
    let args: Vec<String> = env::args().collect();
    let bed_path = &args[1].as_str();
    let ploidy_path = &args[2].as_str();
    let cnv_path = &args[3].as_str();
    let bam_path = args[4].as_str();
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
    let copy_numbers: Vec<usize> = copy_numbers
        .union(&cnv_copy_numbers)
        .filter_map(|x| if *x <= 20 { Some(*x) } else { None })
        .collect();
    eprintln!(
        "Generating compositions for copy numbers {:?}",
        copy_numbers
    );
    let compositions_map = Arc::new(make_compositions_map(&copy_numbers));
    eprintln!("Generated compositions");

    // {
    //     // Even this causes a segfault when reading CRAM...
    //     eprintln!("Entering scope");
    //     let reference = "/Users/maxverbiest/PhD/data/genomes/GRCh38.d1.vd1.fa";
    //     let mut bam = IndexedReader::from_path(bam_path).unwrap();
    //     bam.set_reference(reference).unwrap();
    //     let fetch_request = tr_regions[0].reference_info.get_fetch_definition();
    //     let fetch_result = bam.fetch(fetch_request);
    //     eprintln!("Leaving scope");
    // }
    // eprintln!("Left scope");
    
    eprintln!("Launching {} threads for genotyping", n_threads);
    let chunksize = tr_regions.len() / n_threads + 1;
    tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
        eprintln!("Launched thread {}", current_thread_index().unwrap());
        let htsfile = rhtslib_from_path(bam_path);
        let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
        let c_str = ffi::CString::new(bam_path).unwrap();
        let idx: *mut htslib::hts_idx_t =
            unsafe { htslib::sam_index_load(htsfile, c_str.as_ptr()) };
        if idx.is_null() {
            panic!("Unable to load index for alignment file!");
        }

        let reference = "/Users/maxverbiest/PhD/data/genomes/GRCh38.d1.vd1.fa";
        rhtslib_set_reference(htsfile, reference);

        for tr_region in tr_regions {
            tr_cn_from_cnvs(tr_region, &cnv_regions);

            let fetch_request = tr_region.reference_info.get_fetch_definition_s();
            let itr = rhtslib_fetch_by_str(idx, header, fetch_request.as_bytes());

            fetch_allele_lengths(tr_region, htsfile, itr, flanksize);

            unsafe {
                htslib::hts_itr_destroy(itr);
            }

            estimate_genotype(
                tr_region,
                min_reads_per_allele,
                Arc::clone(&compositions_map),
            );
        }
        unsafe {
            htslib::hts_close(htsfile);
        }
        eprintln!("Finished on thread {}", current_thread_index().unwrap());
    });

    for r in tr_regions {
        println!("{:?}", r);
    }
}


// Functions below that are prefixed with 'rhtslib_' are private functions in
// rust_htslib, and are mostly copied from there to be used in this binary.
// It would be nicer to use rust_htslib's structs (e.g., IndexedReader) for reading 
// alignment files, but those structs caused segfaults when reading CRAM files (not BAM, interestingly)
// when they were dropped, even on trivial tests. -> submit issue on GitHub?
fn rhtslib_from_path<P: AsRef<Path>>(path: P) -> *mut htslib::htsFile {
    let htsfile: *mut htslib::htsFile =
        rhtslib_hts_open(&utils::path_as_bytes(path, true).unwrap(), b"r");
    htsfile
}

fn rhtslib_hts_open(path: &[u8], mode: &[u8]) -> *mut htslib::htsFile {
    let cpath = ffi::CString::new(path).unwrap();
    let c_str = ffi::CString::new(mode).unwrap();
    let ret = unsafe { htslib::hts_open(cpath.as_ptr(), c_str.as_ptr()) };
    if ret.is_null() {
        panic!("Unable to read alignment file!")
    }
    ret
}

fn rhtslib_set_reference<P: AsRef<Path>>(htsfile: *mut htslib::htsFile, path: P) {
    let res = unsafe { bam::set_fai_filename(htsfile, path) };
    if res.is_err() {
        panic!("Problem setting reference file!")
    }
}

fn rhtslib_fetch_by_str(
    idx: *mut htslib::hts_idx_t,
    header: *mut htslib::sam_hdr_t,
    region: &[u8],
) -> *mut htslib::hts_itr_t {
    let rstr = ffi::CString::new(region).unwrap();
    let rptr = rstr.as_ptr();
    let itr = unsafe { htslib::sam_itr_querys(idx, header, rptr) };
    if itr.is_null() {
        panic!("Problem fetching reads from region!")
    }
    itr
}
