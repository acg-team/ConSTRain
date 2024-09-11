//! # Command line interface for `ConSTRain`
use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use log::info;

use crate::utils;

#[derive(Parser)]
#[command(
    name="ConSTRain",
    author,
    version,
    about="Copy number-guided STR allele inference",
    long_about = None,
    propagate_version = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

impl Cli {
    pub fn validate(&self) -> Result<()> {
        match &self.command {
            Commands::Alignment(config) => {
                if let Some(max_norm_depth) = config.max_norm_depth {
                    if config.min_norm_depth > max_norm_depth {
                        bail!("--min-norm-depth cannot be greater than --max-norm-depth")
                    }
                }
                Ok(())
            }
            Commands::VCF(config) => {
                if let Some(max_norm_depth) = config.max_norm_depth {
                    if config.min_norm_depth > max_norm_depth {
                        bail!("--min-norm-depth cannot be greater than --max-norm-depth")
                    }
                }
                Ok(())
            }
        }   
    }

    pub fn get_sample_name(&self) -> Result<String> {
        match &self.command {
            Commands::Alignment(config) => {
                if let Some(name) = &config.sample {
                    Ok(name.clone())
                } else {
                    let name = utils::sample_name_from_path(&config.alignment)?;
                    info!("Inferred sample name from filename: {name}");
                    Ok(name)
                }
            }
            Commands::VCF(config) => Ok(config.sample.clone()),
        }
    }
}

#[derive(Subcommand)]
pub enum Commands {
    /// Extract STR allele lengths from aligment file (SAM/BAM/CRAM) and estimate genotypes
    Alignment(AlignmentArgs),
    /// Extract STR allele lengths from VCF file and estimate genotypes
    VCF(VCFArgs),
}

#[derive(Args)]
pub struct AlignmentArgs {
    /// Input file to extract repeat allele lengths from. Can be SAM/BAM/CRAM.
    #[arg(short, long)]
    pub alignment: String,

    /// File containing chromosome names and their base ploidies. Expected format is JSON
    #[arg(short, long)]
    pub karyotype: String,

    /// File specifying target repeat regions. Expected format is BED3+2
    #[arg(short, long)]
    pub repeats: String,

    /// Copy number variants for this individual. Expected format is BED3+1
    #[arg(long)]
    pub cnvs: Option<String>,    

    /// Size of flanking region around the target repeat that reads need to cover to be considered (both sides)
    #[arg(long, default_value_t = 5)]
    pub flanksize: usize,

    /// Maximum copy number to consider
    #[arg(long, default_value_t = 20)]
    pub max_cn: usize,

    /// Maximum normalised depth of coverage to perform allele length estimation. E.g., `max-norm-depth` of 30. means loci with more than 60 reads for a locus with copy number 2 will be skipped, 90 for copy number 3 etc. (not set by default)
    #[arg(long, value_parser = norm_depth_in_range)]
    pub max_norm_depth: Option<f32>,

    /// Minimum normalised depth of coverage to perform allele length estimation. E.g., `min-norm-depth` of 10. means at least 20 reads are needed for a locus with copy number 2, 30 for copy number 3 etc. (must be at least 1.)
    #[arg(long, default_value_t = 1., value_parser = norm_depth_in_range)]
    pub min_norm_depth: f32,    

    /// Reference genome. Expected format is FASTA (NOT gzipped), index file should exist right next to FASTA. Required if alignment is in CRAM format.
    #[arg(long)]
    pub reference: Option<String>,

    /// Sample name
    #[arg(long)]
    pub sample: Option<String>,

    /// Number of threads to use
    #[arg(long, default_value_t = 1, value_parser = threads_in_range)]
    pub threads: usize,
}

#[derive(Args)]
pub struct VCFArgs {
    /// File containing chromosome names and their base ploidies. Expected format is JSON
    #[arg(short, long)]
    pub karyotype: String,

    /// Sample name
    #[arg(short, long)]
    pub sample: String,

    /// Input file to extract repeat allele lengths from (VCF).
    #[arg(short, long)]
    pub vcf: String,

    /// Copy number variants for this individual. Expected format is BED3+1
    #[arg(long)]
    pub cnvs: Option<String>,

    /// Maximum copy number to consider
    #[arg(long, default_value_t = 20)]
    pub max_cn: usize,

    /// Maximum normalised depth of coverage to perform allele length estimation. E.g., `max-norm-depth` of 30. means loci with more than 60 reads for a locus with copy number 2 will be skipped, 90 for copy number 3 etc. (not set by default)
    #[arg(long, value_parser = norm_depth_in_range)]
    pub max_norm_depth: Option<f32>,    

    /// Minimum normalised depth of coverage to perform allele length estimation. E.g., `min-norm-depth` of 10. means at least 20 reads are needed for a locus with copy number 2, 30 for copy number 3 etc. (must be at least 1.)
    #[arg(long, default_value_t = 1., value_parser = norm_depth_in_range)]
    pub min_norm_depth: f32,

    /// Number of threads to use
    #[arg(long, default_value_t = 1, value_parser = threads_in_range)]
    pub threads: usize,
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

fn norm_depth_in_range(s: &str) -> Result<f32> {
    let dp = s.parse::<f32>().context("Could not parse value passed to --max-norm-depth or --min-norm-depth")?;
    if dp < 1. {
        bail!("--max-norm-depth and --min-norm-depth must be at least 1.")
    }
    Ok(dp)
}
