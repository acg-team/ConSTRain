use std::{collections::HashSet, path::Path, str};

use anyhow::{bail, Context, Result};
use log::info;
use rust_htslib::bcf::{header::Header, Format, Read, Reader, Writer};

use crate::{
    io::RepeatSource,
    karyotype::Karyotype,
    repeat::{RepeatReferenceInfo, TandemRepeat},
    utils::{self, vcf::VcfFilter},
};

pub struct VariantCallFile<P: AsRef<Path>> {
    pub file_path: P,
    pub sample_name: String,
}

impl<P: AsRef<Path>> VariantCallFile<P> {
    pub fn new(file_path: P, sample_name: &str) -> Self {
        Self {
            file_path,
            sample_name: sample_name.into(),
        }
    }
}

impl<P: AsRef<Path>> RepeatSource for VariantCallFile<P> {
    /// Read tandem repeat regions specified in vcf file at `self.file_path` into `repeat_buffer`.
    /// All copy number values that are encoutered while reading repeat regions are added to `copy_number_buffer`.
    fn load_repeats(
        &self,
        karyotype: &Karyotype,
        repeat_buffer: &mut Vec<TandemRepeat>,
        copy_number_buffer: &mut HashSet<usize>,
    ) -> Result<()> {
        let mut bcf = Reader::from_path(&self.file_path).with_context(|| {
            format!(
                "Failed to open VCF file at {}",
                self.file_path.as_ref().display()
            )
        })?;

        let header = bcf.header().to_owned();
        let sample_idx = header
            .sample_id(self.sample_name.as_bytes())
            .with_context(|| {
                format!(
                    "Sample {} not found in file {}",
                    self.sample_name,
                    self.file_path.as_ref().display()
                )
            })?;

        for record in bcf.records() {
            let record = record.with_context(|| {
                format!(
                    "Error reading VCF record in file {}",
                    self.file_path.as_ref().display()
                )
            })?;
            let ref_info = RepeatReferenceInfo::from_bcf_record(&record, &header)?;

            // Generate TandemRepeat using allele lengths from the record, set the genotype to None
            let allele_freqs = utils::vcf::allele_lens_from_record(&record, sample_idx)?;
            let mut tr = TandemRepeat {
                reference_info: ref_info,
                copy_number: 0, // placeholder copy number value
                allele_lengths: allele_freqs,
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

pub struct VariantCallFormatter {
    header: Header,
}

impl VariantCallFormatter {
    pub fn from_vcf_file<P: AsRef<Path>>(vcf: &VariantCallFile<P>) -> Result<Self> {
        let reader = Reader::from_path(&vcf.file_path).with_context(|| {
            format!(
                "Failed to open VCF file at {}",
                vcf.file_path.as_ref().display()
            )
        })?;
        let header = Header::from_template(reader.header());

        Ok(Self { header })
    }
    pub fn from_targets_lengths(
        sample_name: &str,
        targets: &[String],
        lengths: &[u64],
    ) -> Result<Self> {
        if targets.len() == 0 || targets.len() != lengths.len() {
            bail!("Number of targets and number of lengths must be equal and greater than 0");
        }
        let header = utils::vcf::make_bcf_header(targets, lengths, sample_name);

        Ok(Self { header })
    }
    pub fn repeats_to_stdout(&self, tr_regions: &[TandemRepeat]) -> Result<()> {
        let mut vcf = Writer::from_stdout(&self.header, true, Format::Vcf)?;

        for tr_region in tr_regions {
            let record = utils::vcf::repeat_to_bcf_record(&tr_region, &vcf)?;
            vcf.write(&record)?;
        }

        Ok(())
    }
}
