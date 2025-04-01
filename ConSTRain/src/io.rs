use std::collections::{HashMap, HashSet};

use anyhow::Result;
use log::{debug, info};

use crate::{cnv::CopyNumberVariant, karyotype::Karyotype, repeat::TandemRepeat};

pub mod bed;
pub mod json;
pub mod vcf;

pub trait RepeatSource {
    /// Load tandem repeat regions into `repeat_buffer`.
    /// Copy numbers of tandem repeats should be set based on the values specified in `karyotype`. All copy number
    /// values that are observed while reading TR regions should be added to `copy_number_buffer`, which can be used
    /// to keep track of which integers to generate partitions for.
    fn load_repeats(
        &self,
        karyotype: &Karyotype,
        repeat_buffer: &mut Vec<TandemRepeat>,
        copy_number_buffer: &mut HashSet<usize>,
    ) -> Result<()>;
}

pub trait CopyNumberVariantSource {
    /// Read copy number variants into a HashMap where keys are contig
    /// names and values are vectors of CNVs. Implementers should check during the creation of this HashMap
    /// that the CNVs that are being processed (e.g., from a Bed file) are coordinate sorted.
    /// All copy number values that are observed while reading CNVs are added to `copy_number_buffer`,
    /// which can be used to keep track of which integers to generate partitions for.
    fn load_cnvs(
        &self,
        copy_number_buffer: &mut HashSet<usize>,
    ) -> Result<HashMap<String, Vec<CopyNumberVariant>>>;
}

pub fn load_tandem_repeats<T, U>(
    repeat_source: &T,
    karyotype: &str,
    max_cn: usize,
    cnv_source: Option<&U>,
) -> Result<(Vec<TandemRepeat>, Vec<usize>)>
where
    T: RepeatSource,
    U: CopyNumberVariantSource,
{
    let mut repeat_buffer: Vec<TandemRepeat> = Vec::new();
    let mut copy_number_buffer: HashSet<usize> = HashSet::new();
    let karyotype = Karyotype::from_json(karyotype)?;

    repeat_source.load_repeats(&karyotype, &mut repeat_buffer, &mut copy_number_buffer)?;

    if let Some(source) = cnv_source {
        let cnv_map = source.load_cnvs(&mut copy_number_buffer)?;
        for tr_region in &mut repeat_buffer {
            if let Some(cnv_vec) = cnv_map.get(&tr_region.reference_info.seqname) {
                tr_region.set_cn_from_cnvs(cnv_vec);
            };
        }
        info!("Updated TR copy numbers using CNVs");
    }

    let mut observed_copy_numbers: Vec<usize> = copy_number_buffer
        .iter()
        .filter_map(|x| {
            if *x <= max_cn {
                Some(*x)
            } else {
                debug!("Copy number {x} is higher than --max_cn {max_cn}, loci with this CN will not be genotyped");
                None
            }
        })
        .collect();
    observed_copy_numbers.sort();

    Ok((repeat_buffer, observed_copy_numbers))
}
