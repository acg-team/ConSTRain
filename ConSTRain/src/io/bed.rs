use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use anyhow::{bail, Context, Result};
use csv::ReaderBuilder;
use log::info;

use crate::{
    cnv::CopyNumberVariant,
    io::{CopyNumberVariantSource, RepeatSource},
    karyotype::Karyotype,
    repeat::{RepeatReferenceInfo, TandemRepeat},
    utils::vcf::VcfFilter,
};

pub struct BedFile<P: AsRef<Path>> {
    pub file_path: P,
}

impl<P: AsRef<Path>> BedFile<P> {
    pub fn new(file_path: P) -> Self {
        Self { file_path }
    }
}

impl<P: AsRef<Path>> RepeatSource for BedFile<P> {
    /// Read tandem repeat regions from the bed file at `file_path` into `repeat_buffer`.
    /// All copy number values that are encoutered while reading repeat regions are added to `copy_number_buffer`.
    fn load_repeats(
        &self,
        karyotype: &Karyotype,
        repeat_buffer: &mut Vec<TandemRepeat>,
        copy_number_buffer: &mut HashSet<usize>,
    ) -> Result<()> {
        let mut bed_reader = ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(&self.file_path)
            .with_context(|| {
                format!(
                    "Could not read bed file {}",
                    self.file_path.as_ref().display()
                )
            })?;

        for result in bed_reader.deserialize() {
            let ref_info: RepeatReferenceInfo = result.with_context(|| {
                format!(
                    "Failed to deserialize bed record in {}",
                    self.file_path.as_ref().display()
                )
            })?;

            let mut tr = TandemRepeat {
                reference_info: ref_info,
                copy_number: 0, // placeholder copy number value
                allele_lengths: None,
                genotype: None,
                filter: VcfFilter::Pass,
            };

            tr.set_cn_from_karyotpe(karyotype);
            copy_number_buffer.insert(tr.copy_number);
            repeat_buffer.push(tr);
        }

        info!(
            "Read {} TR regions from {}",
            repeat_buffer.len(),
            self.file_path.as_ref().display()
        );
        Ok(())
    }
}

impl<P: AsRef<Path>> CopyNumberVariantSource for BedFile<P> {
    /// Read copy number variants from the bed file at `self.file_path` into a HashMap and return it.
    /// HashMap keys are contig names and values are vectors of CNVs located on that contig.
    /// All copy number values that are encoutered while reader repeat regions are added to `copy_number_buffer`.
    fn load_cnvs(
        &self,
        copy_number_buffer: &mut HashSet<usize>,
    ) -> Result<HashMap<String, Vec<CopyNumberVariant>>> {
        let mut bed_reader = ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(&self.file_path)
            .with_context(|| {
                format!(
                    "Could not read bed file {}",
                    self.file_path.as_ref().display()
                )
            })?;

        let mut cnv_map: HashMap<String, Vec<CopyNumberVariant>> = HashMap::new();
        let mut n = 0;
        for result in bed_reader.deserialize() {
            let cnv: CopyNumberVariant = result.with_context(|| {
                format!(
                    "Failed to deserialize bed record in {}",
                    self.file_path.as_ref().display()
                )
            })?;
            copy_number_buffer.insert(cnv.cn);
            n += 1;
            if let Some(cnv_vec) = cnv_map.get_mut(&cnv.seqname) {
                let prev_cnv = &cnv_vec[cnv_vec.len() - 1];
                if cnv.start < prev_cnv.start {
                    bail!(
                        "CNVs in {} are unsorted. Current: {cnv:?}, previous: {prev_cnv:?}",
                        self.file_path.as_ref().display()
                    );
                } else if cnv.start < prev_cnv.end - 1 {
                    bail!(
                        "Overlapping CNVs in {}. Encountered {cnv:?} and {prev_cnv:?}",
                        self.file_path.as_ref().display()
                    );
                }
                cnv_vec.push(cnv);
            } else {
                cnv_map.insert(cnv.seqname.clone(), vec![cnv]);
            }
        }

        info!("Read {n} CNVS from {}", self.file_path.as_ref().display());

        Ok(cnv_map)
    }
}
