//! # Command line interface for `ConSTRain`
use crate::utils;
use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use log::warn;

#[derive(Parser)]
#[command(
    name="ConSTRain",
    author,
    version,
    about="Copy number-guided STR allele inference",
    long_about = None,
    propagate_version = true
)]
// #[command()]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

impl Cli {
    pub fn get_sample_name(&self) -> Result<String> {
        match &self.command {
            Commands::Alignment {
                repeats: _,
                alignment,
                ploidy: _,
                flanksize: _,
                threads: _,
                reads_per_allele: _,
                sample,
                reference: _,
                cnvs: _,
            } => {
                if let Some(name) = sample {
                    Ok(name.clone())
                } else {
                    let name = utils::sample_name_from_path(alignment)?;
                    warn!("No sample name given. Using inferred sample name: {name}");
                    Ok(name)
                }
            }
            Commands::VCF {} => bail!("VCF subcommand is not implemented"),
        }
    }
}

#[derive(Subcommand)]
pub enum Commands {
    /// Extract allele lengths from aligment file (SAM/BAM/CRAM)
    Alignment {
        /// File specifying target repeat regions. Expected format is BED3+2
        #[arg(short, long)]
        repeats: String,

        /// Input file to extract repeat allele lengths from. Can be SAM/BAM/CRAM.
        #[arg(short, long)]
        alignment: String,

        /// File containing chromosome names and their base ploidies. Expected format is JSON
        #[arg(short, long)]
        ploidy: String,

        /// Size of flanking region around the target repeat that reads need to cover to be considered (both sides)
        #[arg(long, default_value_t = 5)]
        flanksize: usize,

        /// Number of threads to use
        #[arg(long, default_value_t = 1, value_parser = threads_in_range)]
        threads: usize,

        /// Minimum number of reads per copy number to perform allele length estimation. E.g., `reads_per_allele` = 10 means at least 20 reads are needed for a locus with copynumber 2, 30 for copynumber 3 etc.
        #[arg(long, default_value_t = 0)]
        reads_per_allele: usize,

        /// Sample name
        #[arg(long)]
        sample: Option<String>,

        /// Reference genome. Expected format is FASTA, index file should exist right next to FASTA. Required if alignment is in CRAM format.
        #[arg(long)]
        reference: Option<String>,

        /// Copy number variants for this individual. Expected format is BED3+1
        #[arg(long)]
        cnvs: Option<String>,
    },
    /// Extracting allele lengths from VCF is not yet implemented. Please use `ConSTRain alignment`
    VCF {},
}

fn threads_in_range(s: &str) -> Result<usize> {
    let threads = s
        .parse()
        .context("Could not parse value passed to --threads to integer")?;
    if threads < 1 {
        bail!("--threads must be at least 1");
    }
    Ok(threads)
}
