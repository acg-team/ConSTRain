use clap::Parser;
use cn_guided_str_genotying::utils::CopyNumberVariant;
use rayon::{current_thread_index, prelude::*, ThreadPoolBuilder};
use rust_htslib::bam::HeaderView;
use rust_htslib::{bam, htslib, utils};
use std::path::Path;
use std::{ffi, sync::Arc};

use cn_guided_str_genotying::{
    fetch_allele_lengths, genotyping::estimate_genotype, make_partitions_map, tr_cn_from_cnvs,
    utils::io_utils,
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

    /// Minimum number of reads per copy number to perform allele length estimation. E.g., reads_per_allele = 10 means at least 20 reads are needed for a locus with copynumber 2, 30 for copynumber 3 etc.
    #[arg(long, default_value_t = 10)]
    reads_per_allele: usize,
}

fn main() {
    let cli = Cli::parse();

    if cli.threads < 1 {
        panic!("--threads must be at least 1");
    }

    let sample_name = match cli.sample {
        Some(name) => name,
        None => {
            let name = Path::new(&cli.alignment)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap();
            eprintln!(
                "Sample name not specified. Using inferred sample name: {}",
                name
            );
            String::from(name)
        }
    };

    ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    // Currently, io_utils functions return all copy numbers they encounter
    // while reading data. This is then used to create partitions_map.
    // Could instead switch to Arc<RwLock<HashMap<..., ...>>> implementation
    // where partitions_map is updated if a new copy number is encountered
    // in one of the threads.
    let (mut copy_numbers, mut tr_regions) =
        io_utils::trs_from_bed(cli.repeats.as_ref(), cli.ploidy.as_ref());
    eprintln!("Read {} TR regions", tr_regions.len());

    let mut cnv_regions: Vec<CopyNumberVariant> = Vec::new();
    if cli.cnvs.is_some() {
        let (cnv_copy_numbers, mut parsed_cnv_regions) =
            io_utils::cnvs_from_bed(cli.cnvs.unwrap().as_ref());
        cnv_regions.append(&mut parsed_cnv_regions);
        copy_numbers.extend(cnv_copy_numbers.iter());
        eprintln!("Read {} CNVs", cnv_regions.len());
    }
    let copy_numbers: Vec<usize> = copy_numbers
        .iter()
        .filter_map(|x| if *x <= 20 { Some(*x) } else { None })
        .collect();

    eprintln!("Generating partitions for copy numbers {:?}", copy_numbers);
    let partitions_map = Arc::new(make_partitions_map(&copy_numbers));
    eprintln!("Generated partitions");

    eprintln!("Launching {} thread(s) for genotyping", cli.threads);
    let chunksize = tr_regions.len() / cli.threads + 1;
    let alignment_path = cli.alignment.as_str();

    tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
        if cli.threads > 1 {
            eprintln!("Launched thread {}", current_thread_index().unwrap());
        }
        let htsfile = rhtslib_from_path(alignment_path);
        let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
        let c_str = ffi::CString::new(alignment_path).unwrap();
        let idx: *mut htslib::hts_idx_t =
            unsafe { htslib::sam_index_load(htsfile, c_str.as_ptr()) };
        if idx.is_null() {
            panic!("Unable to load index for alignment file!");
        }

        unsafe {
            if htsfile
                .as_ref()
                .expect("Problem accessing htsfile")
                .is_cram()
                != 0
            {
                match &cli.reference {
                    Some(reference) => rhtslib_set_reference(htsfile, reference),
                    None => {
                        panic!("Alignment file is CRAM format but no reference file is specified!")
                    }
                };
            }
        }

        for tr_region in tr_regions {
            if cnv_regions.len() > 0 {
                tr_cn_from_cnvs(tr_region, &cnv_regions);
            }

            let fetch_request = tr_region.reference_info.get_fetch_definition_s();
            let itr = rhtslib_fetch_by_str(idx, header, fetch_request.as_bytes());

            fetch_allele_lengths(tr_region, htsfile, itr, cli.flanksize);

            unsafe {
                htslib::hts_itr_destroy(itr);
            }

            estimate_genotype(tr_region, cli.reads_per_allele, Arc::clone(&partitions_map));
        }
        unsafe {
            // htslib::sam_hdr_destroy(header);
            htslib::hts_close(htsfile);
        }
        if cli.threads > 1 {
            eprintln!("Finished on thread {}", current_thread_index().unwrap());
        }
    });

    // TODO: utility function for getting target names & lengths from alignment file
    let htsfile = rhtslib_from_path(alignment_path);
    let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
    let hview = HeaderView::new(header);

    let mut target_lengths = Vec::<u64>::new();
    for target in hview.target_names().iter() {
        target_lengths.push(hview.target_len(hview.tid(target).unwrap()).unwrap());
    }

    io_utils::trs_to_vcf(
        &tr_regions,
        &hview.target_names(),
        &target_lengths,
        &sample_name,
    );

    unsafe {
        // htslib::sam_hdr_destroy(header);
        htslib::hts_close(htsfile);
    }
}

// Functions below that are prefixed with 'rhtslib_' are private functions in
// rust_htslib and are copied from there to be used in this binary.
// It would be nicer to use rust_htslib's structs (e.g., IndexedReader) for reading
// alignment files, but those structs cause segfaults when reading CRAM files (not BAM, interestingly)
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
