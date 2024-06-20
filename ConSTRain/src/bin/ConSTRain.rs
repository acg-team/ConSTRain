use anyhow::Result;
use clap::Parser;
use constrain::repeat::TandemRepeat;
use constrain::utils::CopyNumberVariant;
use rayon::{prelude::*, ThreadPoolBuilder};
use rust_htslib::{
    bam::{self, HeaderView},
    htslib,
};
use std::collections::HashSet;
use std::{error::Error, ffi, path::Path, sync::Arc};

use constrain::{
    fetch_allele_lengths,
    genotyping::estimate_genotype,
    make_partitions_map, tr_cn_from_cnvs,
    utils::{io_utils, sample_name_from_path},
};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// File specifying target repeat regions. Expected format is BED3+2
    #[arg(long)]
    repeats: String,

    /// File containing chromosome names and their base ploidies. Expected format is JSON
    #[arg(long)]
    ploidy: String,

    /// Sample name
    #[arg(long)]
    sample: Option<String>,

    /// Alignment file to estimate repeat allele lengths from. Can be BAM or CRAM
    #[arg(long)]
    alignment: String,

    /// Copy number variants for the current alignment file. Expected format is BED3+1
    #[arg(long)]
    cnvs: Option<String>,

    /// Reference genome. Expected format is FASTA, index file should exist right next to FASTA. Required if alignment is in CRAM format.
    #[arg(long)]
    reference: Option<String>,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Size of flanking region around the target repeat that reads need to cover to be considered (both sides)
    #[arg(long, default_value_t = 5)]
    flanksize: usize,

    /// Minimum number of reads per copy number to perform allele length estimation. E.g., `reads_per_allele` = 10 means at least 20 reads are needed for a locus with copynumber 2, 30 for copynumber 3 etc.
    #[arg(long, default_value_t = 0)]
    reads_per_allele: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    let sample_name = if let Some(name) = cli.sample {
        name
    } else {
        let name = sample_name_from_path(&cli.alignment)?;
        eprintln!("Sample name not specified. Using inferred sample name: {name}");
        name
    };

    // Currently, io_utils functions return all copy numbers they encounter
    // while reading data. This is then used to create partitions_map.
    // Could instead switch to Arc<RwLock<HashMap<..., ...>>> implementation
    // where partitions_map is updated if a new copy number is encountered
    // in one of the threads.
    let mut copy_numbers: HashSet<usize> = HashSet::new();
    let mut tr_regions: Vec<TandemRepeat> = Vec::new();
    io_utils::trs_from_bed(
        &cli.repeats,
        &cli.ploidy,
        &mut tr_regions,
        &mut copy_numbers,
    )?;
    eprintln!("Read {} TR regions", tr_regions.len());

    let mut cnv_regions: Vec<CopyNumberVariant> = Vec::new();
    if let Some(cnv_file) = cli.cnvs {
        io_utils::cnvs_from_bed(&cnv_file, &mut cnv_regions, &mut copy_numbers)?;
        eprintln!("Read {} CNVs", cnv_regions.len());
    }

    let copy_numbers: Vec<usize> = copy_numbers
        .iter()
        .filter_map(|x| if *x <= 20 { Some(*x) } else { None }) // ConSTRain supports copy numbers up to 20 for now
        .collect();

    eprintln!("Generating partitions for copy numbers {copy_numbers:?}");
    let partitions_map = Arc::new(make_partitions_map(&copy_numbers));
    eprintln!("Generated partitions");

    assert!(cli.threads >= 1, "--threads must be at least 1");
    eprintln!("Launching {} thread(s) for genotyping", cli.threads);
    ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()?;
    let chunksize = tr_regions.len() / cli.threads + 1;
    tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
        // Main work happens in this parallel iterator
        //
        // For each thread, we instantiate a separate `htsfile`, `header`, and `idx`.
        // If anything goes wrong in this process we panic, since these are absolutely essential.
        // After the setup, we iterate over tandem repeat regions, create an iterator for each region,
        // fetch allele lengths from the mapped reads, and estimate the underlying genotype from the 
        // alle length distribution. If something goes wrong here, we just log the error and continue 
        // to the next tandem repeat region.
        let tidx = rayon::current_thread_index().unwrap_or(0);
        if cli.threads > 1 {
            eprintln!("Launched thread {tidx}");
        }

        let htsfile = rhtslib_from_path(&cli.alignment);
        let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
        // we panic if we can't convert the alignment file name to a C string
        let c_str = ffi::CString::new(cli.alignment.as_str()).expect("Internal 0 byte contained in alignment file name");
        let idx: *mut htslib::hts_idx_t =
            unsafe { htslib::sam_index_load(htsfile, c_str.as_ptr()) };
        assert!(!idx.is_null(), "Unable to load index for alignment file!");

        unsafe {
            let is_cram = htsfile
                .as_ref()
                .expect("Problem accessing htsfile")
                .is_cram()
                != 0;
            if is_cram {
                match &cli.reference {
                    Some(reference) => rhtslib_set_reference(htsfile, reference),
                    None => {
                        panic!("Alignment file is CRAM format but no reference file is specified!")
                    }
                };
            }
        }

        for tr_region in tr_regions {
            let fetch_request = tr_region.reference_info.get_fetch_definition_s();
            if !cnv_regions.is_empty() {
                if let Err(e) = tr_cn_from_cnvs(tr_region, &cnv_regions) {
                    eprintln!("Thread {tidx}: Error setting copy number for region {fetch_request}: {e:?}");
                    continue
                };
            }

            let itr = rhtslib_fetch_by_str(idx, header, fetch_request.as_bytes());
            if let Err(e) = fetch_allele_lengths(tr_region, htsfile, itr, cli.flanksize) {
                eprintln!("Thread {tidx}: Error fetching reads for region {fetch_request}: {e:?}");
                // destroy iterator and continue to the next repeat region
                unsafe {
                    htslib::hts_itr_destroy(itr);
                }
                continue
            };
            // destroy iterator
            unsafe {
                htslib::hts_itr_destroy(itr);
            }

            if let Err(e) =
                estimate_genotype(tr_region, cli.reads_per_allele, Arc::clone(&partitions_map))
            {
                eprintln!("Thread {tidx}: Could not estimate genotype for repeat {fetch_request}: {e:?}");
                continue
            }
        }
        unsafe {
            htslib::hts_close(htsfile);
        }
        if cli.threads > 1 {
            eprintln!("Finished on thread {tidx}");
        }
    });

    let (target_names, target_lengths) = get_target_names_lengths(cli.alignment.as_str())?;
    io_utils::trs_to_vcf(&tr_regions, &target_names, &target_lengths, &sample_name)?;

    Ok(())
}

fn get_target_names_lengths(
    alignment_path: &str,
) -> Result<(Vec<String>, Vec<u64>), Box<dyn Error>> {
    let htsfile = rhtslib_from_path(alignment_path);
    let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
    let hview = HeaderView::new(header);
    unsafe {
        htslib::hts_close(htsfile);
    }

    let mut target_names = Vec::<String>::new();
    let mut target_lengths = Vec::<u64>::new(); // tlens are u64 in rust_htslib

    for target in &hview.target_names() {
        let tid = hview
            .tid(target)
            .ok_or("Could not get target ID from header")?;
        let tlen = hview
            .target_len(tid)
            .ok_or("Could not get target length from header")?;
        let tname = std::str::from_utf8(target)?.to_owned();

        target_lengths.push(tlen);
        target_names.push(tname);
    }

    Ok((target_names, target_lengths))
}

// Functions below that are prefixed with 'rhtslib_' are private functions in
// rust_htslib and are copied from there to be used in this binary.
// It would be nicer to use rust_htslib's structs (e.g., IndexedReader) for reading
// alignment files, but those structs cause segfaults when reading CRAM files (not BAM, interestingly)
// when they were dropped, even on trivial tests. -> submit issue on GitHub?
fn rhtslib_from_path<P: AsRef<Path>>(path: P) -> *mut htslib::htsFile {
    let htsfile: *mut htslib::htsFile = rhtslib_hts_open(
        &rust_htslib::utils::path_as_bytes(path, true).unwrap(),
        b"r",
    );
    htsfile
}

fn rhtslib_hts_open(path: &[u8], mode: &[u8]) -> *mut htslib::htsFile {
    let cpath = ffi::CString::new(path).unwrap();
    let c_str = ffi::CString::new(mode).unwrap();
    let ret = unsafe { htslib::hts_open(cpath.as_ptr(), c_str.as_ptr()) };
    assert!(!ret.is_null(), "Unable to read alignment file!");

    ret
}

fn rhtslib_set_reference<P: AsRef<Path>>(htsfile: *mut htslib::htsFile, path: P) {
    let res = unsafe { bam::set_fai_filename(htsfile, path) };
    assert!(res.is_ok(), "Problem setting reference file!");
}

fn rhtslib_fetch_by_str(
    idx: *mut htslib::hts_idx_t,
    header: *mut htslib::sam_hdr_t,
    region: &[u8],
) -> *mut htslib::hts_itr_t {
    let rstr = ffi::CString::new(region).unwrap();
    let rptr = rstr.as_ptr();
    let itr = unsafe { htslib::sam_itr_querys(idx, header, rptr) };
    assert!(!itr.is_null(), "Problem fetching reads from region!");

    itr
}
