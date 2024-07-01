use anyhow::Result;
use clap::Parser;
use constrain::{self, cli::Cli, utils::io_utils};
use env_logger::{Builder, Env};
use log::info;
use rayon::{prelude::*, ThreadPoolBuilder};
use std::sync::Arc;

fn main() -> Result<()> {
    // Initialize the logger. If the log level is not set via `RUST_LOG`, set it to 'info' by default
    Builder::from_env(Env::default().default_filter_or("info")).init();

    // parse command line and validate inputs where possible
    let config = Cli::parse();
    let sample_name = config.get_sample_name()?;

    // parse tandem repeats from bed file and store observed copy numbers
    let (mut tr_regions, observed_copy_numbers) =
        io_utils::parse_tandem_repeats(&config.repeats, &config.ploidy, config.cnvs.as_deref())?;
    
    // generate partitions relevant to genotyping tandem repeats with observed copy numbers
    info!("Generating partitions for copy numbers {observed_copy_numbers:?}");
    let partitions_map = Arc::new(constrain::make_partitions_map(&observed_copy_numbers));

    ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()?;
    let chunksize = tr_regions.len() / config.threads + 1;

    info!("Starting genotyping");
    tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
        // Main work happens in this parallel iterator
        let tidx = rayon::current_thread_index().unwrap_or(0);
        constrain::run(
            tr_regions,
            &partitions_map,
            &config.alignment,
            config.reference.as_deref(),
            config.flanksize,
            config.reads_per_allele,
            tidx,
        );
    });
    info!("Finished genotyping");

    // get contig names and lengths, write variant calls to stdout
    let (target_names, target_lengths) = io_utils::tnames_tlens_from_header(&config.alignment)?;
    io_utils::trs_to_vcf(&tr_regions, &target_names, &target_lengths, &sample_name)?;

    Ok(())
}
