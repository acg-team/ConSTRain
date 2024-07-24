use std::collections::HashSet;

use anyhow::{bail, Result};
use log::{debug, info};

use crate::{karyotype::Karyotype, repeat::TandemRepeat};

pub mod bed;
pub mod json;
pub mod vcf;

#[derive(Debug)]
pub enum InputFileType {
    Alignment(String),
    VCF(String),
}

pub fn parse_tandem_repeats(
    repeats: InputFileType,
    ploidy: &str,
    max_cn: usize,
    cnvs: Option<&str>,
    sample: Option<&str>,
) -> Result<(Vec<TandemRepeat>, Vec<usize>)> {
    let mut observed_copy_numbers: HashSet<usize> = HashSet::new();
    let mut tr_regions: Vec<TandemRepeat> = Vec::new();
    let karyotype = Karyotype::from_json(&ploidy)?;

    #[allow(unreachable_patterns)]
    match repeats {
        InputFileType::Alignment(path) => {
            bed::read_trs(
                &path,
                &karyotype,
                &mut tr_regions,
                &mut observed_copy_numbers,
            )?;
        }
        InputFileType::VCF(path) => {
            if let Some(sample) = sample {
                vcf::read_trs(
                    &path,
                    &karyotype,
                    &mut tr_regions,
                    &mut observed_copy_numbers,
                    &sample,
                )?;
            } else {
                bail!("Sample name needs to be set to read repeats from VCF input");
            }
        }
        _ => bail!("Parsing repeats from filetype {repeats:?} is not supported"),
    }

    if let Some(cnv_file) = cnvs {
        let cnv_map = bed::read_cnvs(cnv_file, &mut observed_copy_numbers)?;
        for tr_region in &mut tr_regions {
            if let Some(cnv_vec) = cnv_map.get(&tr_region.reference_info.seqname) {
                if let Err(_) = tr_region.set_cn_from_cnvs(cnv_vec) {
                    debug!(
                        "Could not set copy number of {} from CNV, skipping locus",
                        tr_region.reference_info.get_fetch_definition_s()
                    );
                    tr_region.skip = true;
                    continue;
                };
            };
        }
        info!("Updated TR copy numbers using CNVs");
    }

    let mut observed_copy_numbers: Vec<usize> = observed_copy_numbers
        .iter()
        .filter_map(|x| {
            if *x <= max_cn {
                Some(*x)
            } else {
                debug!("Copy number {x} is outside of the currently supported range of 1-{max_cn}, loci with this CN will not be genotyped");
                None
            }
        })
        .collect();
    observed_copy_numbers.sort();

    Ok((tr_regions, observed_copy_numbers))
}
