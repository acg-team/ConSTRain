use anyhow::{Context, Result};
use clap::Parser;
use constrain::repeat::TandemRepeat;
use constrain::utils::CopyNumberVariant;
use rayon::{current_thread_index, prelude::*, ThreadPoolBuilder};
use rust_htslib::{bam::HeaderView, htslib};
use std::collections::HashSet;
use std::{error::Error, sync::Arc};

use constrain::{
    self, rhtslib_reimplement,
    utils::{self, io_utils},
};

#[derive(Parser)]
#[command(name = "ConSTRain")]
#[command(version)]
#[command(about = "Copy number-guided STR allele inference", long_about = None)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// File specifying target repeat regions. Expected format is BED3+2
    #[arg(short, long)]
    repeats: String,

    /// File containing chromosome names and their base ploidies. Expected format is JSON
    #[arg(short, long)]
    ploidy: String,

    /// Alignment file to estimate repeat allele lengths from. Can be BAM or CRAM
    #[arg(short, long)]
    alignment: String,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    /// Size of flanking region around the target repeat that reads need to cover to be considered (both sides)
    #[arg(short, long, default_value_t = 5)]
    flanksize: usize,

    /// Minimum number of reads per copy number to perform allele length estimation. E.g., `reads_per_allele` = 10 means at least 20 reads are needed for a locus with copynumber 2, 30 for copynumber 3 etc.
    #[arg(long, default_value_t = 0)]
    reads_per_allele: usize,

    /// Sample name
    #[arg(long)]
    sample: Option<String>,

    /// Copy number variants for this individual. Expected format is BED3+1
    #[arg(long)]
    cnvs: Option<String>,

    /// Reference genome. Expected format is FASTA, index file should exist right next to FASTA. Required if alignment is in CRAM format.
    #[arg(long)]
    reference: Option<String>,
}

fn main() -> Result<(), Box<dyn Error>> {
    // parse command line and validate inputs where possible
    let config = Cli::parse();

    let sample_name = if let Some(name) = config.sample {
        name
    } else {
        let name = utils::sample_name_from_path(&config.alignment)?;
        eprintln!("Sample name not specified. Using inferred sample name: {name}");
        name
    };
    assert!(config.threads >= 1, "--threads must be at least 1");

    let (mut tr_regions, cnv_regions, observed_copy_numbers) =
        read_tandem_repeats(&config.repeats, &config.ploidy, config.cnvs.as_ref())?;
    eprintln!("Generating partitions for copy numbers {observed_copy_numbers:?}");
    let partitions_map = Arc::new(constrain::make_partitions_map(&observed_copy_numbers));

    ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()?;
    let chunksize = tr_regions.len() / config.threads + 1;

    tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
        // Main work happens in this parallel iterator
        let tidx = current_thread_index().unwrap_or(0);
        constrain::run(
            tr_regions,
            &cnv_regions,
            &partitions_map,
            &config.alignment,
            config.reference.as_ref(),
            config.flanksize,
            config.reads_per_allele,
            tidx,
        );
    });

    let (target_names, target_lengths) = get_target_names_lengths(&config.alignment)?;
    io_utils::trs_to_vcf(&tr_regions, &target_names, &target_lengths, &sample_name)?;

    Ok(())
}

fn get_target_names_lengths(alignment_path: &str) -> Result<(Vec<String>, Vec<u64>)> {
    let htsfile = rhtslib_reimplement::rhtslib_from_path(alignment_path)?;
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
            .context("Could not get target ID from header")?;
        let tlen = hview
            .target_len(tid)
            .context("Could not get target length from header")?;
        let tname = std::str::from_utf8(target)?.to_owned();

        target_lengths.push(tlen);
        target_names.push(tname);
    }

    Ok((target_names, target_lengths))
}

fn read_tandem_repeats(
    repeats: &str,
    ploidy: &str,
    cnvs: Option<&String>,
) -> Result<(Vec<TandemRepeat>, Vec<CopyNumberVariant>, Vec<usize>)> {
    // Currently, io_utils functions return all copy numbers they encounter
    // while reading data. This is then used to create partitions_map.
    // Could instead switch to Arc<RwLock<HashMap<..., ...>>> implementation
    // where partitions_map is updated if a new copy number is encountered
    // in one of the threads.
    let mut observed_copy_numbers: HashSet<usize> = HashSet::new();
    let mut tr_regions: Vec<TandemRepeat> = Vec::new();
    io_utils::trs_from_bed(repeats, ploidy, &mut tr_regions, &mut observed_copy_numbers)?;
    eprintln!("Read {} TR regions", tr_regions.len());

    let mut cnv_regions: Vec<CopyNumberVariant> = Vec::new();
    if let Some(cnv_file) = cnvs {
        io_utils::cnvs_from_bed(&cnv_file, &mut cnv_regions, &mut observed_copy_numbers)?;
        eprintln!("Read {} CNVs", cnv_regions.len());
    }

    let mut observed_copy_numbers: Vec<usize> = observed_copy_numbers
        .iter()
        .filter_map(|x| if 0 < *x && *x <= 20 { Some(*x) } else { None }) // ConSTRain supports copy numbers from 1 to 20 for now
        .collect();
    observed_copy_numbers.sort();

    Ok((tr_regions, cnv_regions, observed_copy_numbers))
}
