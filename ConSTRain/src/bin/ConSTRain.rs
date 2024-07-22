use anyhow::Result;
use clap::Parser;
use constrain::{
    self,
    cli::{Cli, Commands},
    genotyping, io, utils,
};
use env_logger::{Builder, Env};
use log::{debug, info};
use rayon::{prelude::*, ThreadPoolBuilder};
use std::sync::Arc;

fn main() -> Result<()> {
    // Initialize the logger. If the log level is not set via `RUST_LOG`, default to 'info'
    Builder::from_env(Env::default().default_filter_or("info")).init();

    // parse command line and validate inputs where possible
    let config = Cli::parse();
    let sample_name = config.get_sample_name()?;
    debug!("Parsed command line");

    // parse tandem repeats from bed file and store observed copy numbers
    match config.command {
        Commands::Alignment(args) => {
            // read tandem repeats from bedfile. Copy numbers will be set based on `ploidy`, and
            // optionally updated from `cnvs`, if it was provided
            let (mut tr_regions, observed_copy_numbers) =
                io::bed::parse_tandem_repeats(&args.repeats, &args.ploidy, args.cnvs.as_deref())?;

            // generate partitions relevant to genotyping tandem repeats with observed copy numbers
            let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

            debug!("Spawning {} threads", args.threads);
            ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build_global()?;
            let chunksize = tr_regions.len() / args.threads + 1;

            info!("Starting genotyping");
            tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
                // Main work happens in this parallel iterator
                let tidx = rayon::current_thread_index().unwrap_or(0);
                constrain::run(
                    tr_regions,
                    &partitions_map,
                    &args.alignment,
                    args.reference.as_deref(),
                    args.flanksize,
                    args.reads_per_allele,
                    tidx,
                )
                .expect("Error during genotyping");
            });
            info!("Finished genotyping");

            // get contig names and lengths, write variant calls to stdout
            let (target_names, target_lengths) = utils::tnames_tlens_from_header(&args.alignment)?;
            io::vcf::write(&tr_regions, &target_names, &target_lengths, &sample_name)?;
        }
        Commands::VCF(args) => {
            panic!("VCF subcommand is not implemented");
        }
    };

    Ok(())
}
