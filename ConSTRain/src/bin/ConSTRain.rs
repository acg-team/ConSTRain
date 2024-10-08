use anyhow::Result;
use clap::Parser;
use constrain::{
    self,
    cli::{Cli, Commands},
    genotyping,
    io::{self, InputFileType},
    utils,
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
    config.validate()?;
    let sample_name = config.get_sample_name()?;
    debug!("Parsed command line");

    match config.command {
        Commands::Alignment(args) => {
            // read tandem repeats from bedfile. Copy numbers will be set based on `karyotype` and
            // optionally updated from `cnvs` if it was provided
            let (mut tr_regions, observed_copy_numbers) = io::load_tandem_repeats(
                InputFileType::Alignment(args.repeats.to_string()),
                &args.karyotype,
                args.max_cn,
                args.cnvs.as_deref(),
                args.sample.as_deref(),
            )?;

            // generate partitions relevant to genotyping tandem repeats with observed copy numbers
            let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

            debug!("Spawning {} thread(s) for genotyping", args.threads);
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
                    args.min_norm_depth,
                    args.max_norm_depth,
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
            let (mut tr_regions, observed_copy_numbers) = io::load_tandem_repeats(
                InputFileType::VCF(args.vcf.to_string()),
                &args.karyotype,
                args.max_cn,
                args.cnvs.as_deref(),
                Some(&args.sample),
            )?;

            // generate partitions relevant to genotyping tandem repeats with observed copy numbers
            let partitions_map = Arc::new(genotyping::make_partitions_map(&observed_copy_numbers));

            debug!("Spawning {} thread(s) for genotyping", args.threads);
            ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build_global()?;
            let chunksize = tr_regions.len() / args.threads + 1;

            info!("Starting genotyping");
            tr_regions.par_chunks_mut(chunksize).for_each(|tr_regions| {
                // Main work happens in this parallel iterator
                let tidx = rayon::current_thread_index().unwrap_or(0);
                constrain::run_vcf(tr_regions, &partitions_map, args.min_norm_depth, args.max_norm_depth, tidx)
                    .expect("Error during genotyping");
            });

            info!("Finished genotyping");

            io::vcf::write_reuse_header(&tr_regions, &args.vcf)?;
        }
    };

    Ok(())
}
