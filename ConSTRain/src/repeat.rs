//! # Structs to represent genomic tandem repeats
//!
//! Module containing structs to represent tandem repeats.
//! [`TandemRepeat`] is the struct that represents an observation of a tandem
//! repeat in an alignment, included observed allele lengths across reads and the
//! copy number in a specific sample. It is also associated with a [`RepeatReferenceInfo`] struct,
//! which encodes how a repetetive region looks in the reference genome.
//! `ConSTRain` assumes that TRs are perfect tandem repetitions of a given nucleotide motif.
use std::{collections::HashMap, str};

use anyhow::{Context, Result};
use log::debug;
use ndarray::prelude::*;
use rust_htslib::bcf::{header::HeaderView, Record};

use crate::{
    cnv::CopyNumberVariant,
    karyotype::Karyotype,
    utils::{self, vcf::VcfFilter},
};

/// `TandemRepeat` represents an observation of a tandem in an alignment.
/// The representation of the locus in the reference genome is contained in `reference_info`. The
/// copy number of the tandem repeat in the current sample is given by `copy_number`. `allele_lengths`
/// contains the observed allele lengths for the repeat in the current alignment. `genotype` is the underlying
/// genotype that is inferred to have given rise to the observed allele length distribution.
#[derive(Debug)]
pub struct TandemRepeat {
    pub reference_info: RepeatReferenceInfo,
    pub copy_number: usize,
    pub allele_lengths: Option<HashMap<i64, f32>>,
    pub genotype: Option<Vec<(i64, f32)>>,
    // pub skip: bool, // if false, do not estimate genotype for this TR
    pub filter: VcfFilter,
}

impl TandemRepeat {
    pub fn has_coverage(&self) -> bool {
        self.allele_lengths.is_some()
    }
    pub fn set_allele_lengths(&mut self, new_allele_lengths: Option<HashMap<i64, f32>>) {
        self.allele_lengths = new_allele_lengths
    }
    pub fn set_cn(&mut self, new_cn: usize) {
        self.copy_number = new_cn;
    }
    /// Search for a CNV in `cnv_regions` that overlaps with this tandem repeat.
    /// If there is such an overlap: update `self.copy_number`.
    /// Assumes that all entries in `cnv_regions` are on the same contig as this tandem repeat.
    /// Furthermore, we assume that `cnv_regions` is coordinate sorted. **These assumptions
    /// are not checked here. Improper input may result in the wrong copy number being set!**
    pub fn set_cn_from_cnvs(&mut self, cnv_regions: &[CopyNumberVariant]) {
        let region_len = self.reference_info.end - self.reference_info.start;
        for cnv in cnv_regions {
            if cnv.start > self.reference_info.end {
                // The current CNV region lies beyond the STR on the contig. Since CNVs are coordinate
                // sorted, we can be confident that no other CNV overlaps the TR.
                break;
            }

            let Ok(overlap) = utils::range_overlap(
                self.reference_info.start,
                self.reference_info.end - 1,
                cnv.start,
                cnv.end - 1,
            ) else {
                debug!(
                    "Error calculating overlap between {} and {:?}",
                    self.reference_info.get_fetch_definition_s(),
                    cnv
                );
                continue;
            };

            if overlap == region_len {
                debug!(
                    "Updating CN from {} to {} for {}",
                    self.copy_number,
                    cnv.cn,
                    self.reference_info.get_fetch_definition_s()
                );
                self.set_cn(cnv.cn);
                // TRs can intersect with at most one CNV, we found a hit so we can return
                break;
            } else if overlap > 0 {
                // TR partially overlaps CNA, impossible to set sensible CN for TR, so we update the filter type
                debug!(
                    "Locus {} only partially overlaps a CNA and will not be genotyped",
                    self.reference_info.get_fetch_definition_s()
                );
                self.filter = VcfFilter::CnMissing;
                break;
            }
        }
    }
    pub fn set_cn_from_karyotpe(&mut self, karyotype: &Karyotype) {
        let cn = karyotype.get_ploidy(&self.reference_info.seqname);
        if let Some(cn) = cn {
            self.set_cn(cn);
        } else {
            debug!(
                "Could not set copy number of {} from karyotype, skipping locus",
                self.reference_info.get_fetch_definition_s()
            );
            self.filter = VcfFilter::CnMissing;
        }
    }
    pub fn allele_freqs_as_tuples(&self) -> Vec<(i64, f32)> {
        match &self.allele_lengths {
            Some(allele_lengths) => allele_lengths.iter().map(|i| (*i.0, *i.1)).collect(),
            None => Vec::new(),
        }
    }
    pub fn allele_freqs_as_ndarrays(
        &self,
        sort_by: Option<&str>,
    ) -> (Array<i64, Dim<[usize; 1]>>, Array<f32, Dim<[usize; 1]>>) {
        let mut count_vec = self.allele_freqs_as_tuples();

        if let Some(by) = sort_by {
            if by == "freq" {
                count_vec.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            } else if by == "len" {
                count_vec.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            } else {
                eprintln!("If provided, sort_by must be 'freq' or 'len'. Skipping sort");
            }
        }

        let allele_lengths = Array::<i64, _>::from_vec(count_vec.iter().map(|a| a.0).collect());
        let counts = Array::<f32, _>::from_vec(count_vec.iter().map(|a| a.1).collect());

        (allele_lengths, counts)
    }
    pub fn gt_as_allele_lengths(&self) -> Vec<i64> {
        match &self.genotype {
            Some(genotype) => {
                let mut allele_lens: Vec<i64> = Vec::new();
                for allele in genotype {
                    for _ in 0..allele.1 as usize {
                        allele_lens.push(allele.0);
                    }
                }
                allele_lens
            }
            None => Vec::new(),
        }
    }
    pub fn get_n_mapped_reads(&self) -> Option<usize> {
        self.allele_lengths
            .as_ref()
            .map(|allele_lengths| allele_lengths.values().sum::<f32>() as usize)
    }

    pub fn get_norm_depth(&self) -> Option<f32> {
        self.get_n_mapped_reads()
            .map(|dp| dp as f32 / self.copy_number as f32)
    }
}

/// `RepeatReferenceInfo` stores information encoding how a [`TandemRepeat`] is represented
/// in the reference genome. The `start` and `end` are encoded in the BED coordinate system:
/// 0-based, half-open. `unit` represents the repeating nucleotide motif, and `period` the length
/// of `unit` (convenient to store separately so we don't need to do `unit.len()` all the time).
#[derive(Debug, serde::Deserialize)]
pub struct RepeatReferenceInfo {
    pub seqname: String,
    // Tandem repeat struct follows 0-based half-open coordinate system: [start, end)
    pub start: i64,
    pub end: i64,
    pub period: i64,
    pub unit: String,
}

impl RepeatReferenceInfo {
    pub fn from_bcf_record(record: &Record, header: &HeaderView) -> Result<Self> {
        let rid = record.rid().context("Failed to get record rid")?;
        let contig = str::from_utf8(header.rid2name(rid)?)?.to_string();
        let end = utils::vcf::get_info_int(record, "END")?;
        let period = utils::vcf::get_info_int(record, "PERIOD")?;
        let unit = utils::vcf::get_info_str(record, "RU")?;

        Ok(RepeatReferenceInfo {
            seqname: contig,
            start: record.pos(),
            end,
            period,
            unit,
        })
    }
    pub fn get_fetch_definition(&self) -> (&str, i64, i64) {
        (self.seqname.as_str(), self.start, self.end)
    }
    pub fn get_fetch_definition_s(&self) -> String {
        format!("{}:{}-{}", self.seqname, self.start, self.end)
    }
    pub fn get_reference_len(&self) -> i64 {
        assert!(
            (self.end - self.start) % self.period == 0,
            "Start and end positions for repeat '{}' are not compatible with the period!",
            self.get_fetch_definition_s()
        );

        (self.end - self.start) / self.period
    }
}
