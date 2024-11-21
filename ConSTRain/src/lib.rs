//! # ConSTRain
//!
//! This library serves as the backbone for the [ConSTRain binary](https://github.com/acg-team/ConSTRain),
//! which is developed and maintained by the Applied Computational Genomics Team of the
//! [Bioinformatics Centre](https://www.zhaw.ch/en/lsfm/institutes-centres/icls/bioinformatics/)
//! at the Zürich University of Applied Sciences, Switzerland.
pub mod cli;
pub mod cnv;
pub mod genotyping;
pub mod io;
pub mod karyotype;
pub mod repeat;
pub mod rhtslib_reimplements;
pub mod utils;

use anyhow::{Context, Result};
use genotyping::PartitionMap;
use log::{debug, trace};
use rust_htslib::{
    bam::{ext::BamRecordExtensions, record::CigarStringView, Record},
    htslib::{self, htsFile},
};
use std::{collections::HashMap, ffi, sync::Arc};
use utils::VcfFilter;

use crate::{repeat::TandemRepeat, utils::cigar};

/// The main work of ConSTRain happens in this `run` function.
/// It is meant to be called from inside a rayon parallel iterator.
/// For each thread, we instantiate a separate `htsfile`, `header`, and `idx`.
/// If anything goes wrong in this process we panic, since these are absolutely essential.
/// After the setup, we iterate over tandem repeat regions, create an iterator for each region,
/// fetch allele lengths from the mapped reads, and estimate the underlying genotype from the
/// alle length distribution. If something goes wrong here, we just log the error and continue
/// to the next tandem repeat region.
pub fn run(
    tr_regions: &mut [TandemRepeat],
    partitions_map: &Arc<PartitionMap>,
    alignment: &str,
    reference: Option<&str>,
    flanksize: usize,
    min_norm_depth: f32,
    max_norm_depth: Option<f32>,
    tidx: usize,
) -> Result<()> {
    trace!("Launching thread {tidx}");

    let (htsfile, idx, header) = thread_setup(alignment, reference)
        .with_context(|| format!("Error during setup on thread {tidx}"))?;

    for tr_region in tr_regions {
        if !matches!(tr_region.filter, VcfFilter::Pass) {
            continue;
        }

        let fetch_request = tr_region.reference_info.get_fetch_definition_s();

        let Ok(itr) =
            rhtslib_reimplements::rhtslib_fetch_by_str(idx, header, fetch_request.as_bytes())
        else {
            debug!("Error fetching reads, skipping locus {fetch_request}");
            continue;
        };

        if let Err(e) = extract_allele_lengths(tr_region, htsfile, itr, flanksize) {
            debug!("Error extracting allele lengths, skipping locus {fetch_request}: {e:?}");
            tr_region.filter = VcfFilter::Undef;
            // destroy iterator and continue to the next repeat region
            unsafe {
                htslib::hts_itr_destroy(itr);
            }
            continue;
        };
        // destroy iterator
        unsafe {
            htslib::hts_itr_destroy(itr);
        }

        let pmap = Arc::clone(partitions_map); 
        if let Err(e) = genotyping::tr_region_precheck(tr_region, min_norm_depth, max_norm_depth, &pmap) {
            debug!(
                "Locus {} failed precheck: {e:?}",
                tr_region.reference_info.get_fetch_definition_s()
            );
            continue;
        }

        if let Err(e) = genotyping::estimate_genotype(
            tr_region,
            pmap.get(&tr_region.copy_number).unwrap(), // unwrap here because we check in `tr_region_precheck()` whether entry exists
        ) {
            debug!("Could not estimate genotype for locus {fetch_request}: {e:?}");
            continue;
        }
    }
    unsafe {
        htslib::hts_close(htsfile);
    }

    trace!("Finished on thread {tidx}");
    Ok(())
}

pub fn run_vcf(
    tr_regions: &mut [TandemRepeat],
    partitions_map: &Arc<PartitionMap>,
    min_norm_depth: f32,
    max_norm_depth: Option<f32>,
    tidx: usize,
) -> Result<()> {
    trace!("Launching thread {tidx}");

    for tr_region in tr_regions {
        if !matches!(tr_region.filter, VcfFilter::Pass) {
            continue;
        }
        let pmap = Arc::clone(partitions_map);
        if let Err(e) = genotyping::tr_region_precheck(tr_region, min_norm_depth, max_norm_depth, &pmap) {
            debug!(
                "Locus {} failed precheck: {e:?}",
                tr_region.reference_info.get_fetch_definition_s()
            );
            continue;
        }
        if let Err(e) = genotyping::estimate_genotype(
            tr_region,
            pmap.get(&tr_region.copy_number).unwrap(), // unwrap here because we check in `tr_region_precheck()` whether entry exists
        ) {
            debug!(
                "Could not estimate genotype for locus {}: {e:?}",
                tr_region.reference_info.get_fetch_definition_s()
            );
            continue;
        }
    }

    trace!("Finished on thread {tidx}");
    Ok(())
}

fn thread_setup(
    alignment_path: &str,
    reference: Option<&str>,
) -> Result<(*mut htsFile, *mut htslib::hts_idx_t, *mut htslib::sam_hdr_t)> {
    let htsfile = rhtslib_reimplements::rhtslib_from_path(alignment_path)?;
    let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
    let c_str = ffi::CString::new(alignment_path)
        .context("Encountered internal 0 byte in alignment file name")?;
    let idx: *mut htslib::hts_idx_t = unsafe { htslib::sam_index_load(htsfile, c_str.as_ptr()) };
    assert!(!idx.is_null(), "Unable to load index for alignment file!");

    unsafe {
        let is_cram = htsfile
            .as_ref()
            .with_context(|| "Problem acessing htsfile")?
            .is_cram()
            != 0;

        if is_cram {
            let reference = reference
                .context("Alignment file is CRAM format but no reference file is specified!")?;
            rhtslib_reimplements::rhtslib_set_reference(htsfile, reference)?;
        }
    }

    Ok((htsfile, idx, header))
}

pub fn extract_allele_lengths(
    tr_region: &mut TandemRepeat,
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    flank: usize,
) -> Result<()> {
    let mut allele_lengths: HashMap<i64, f32> = HashMap::new();
    let mut record = Record::new();
    while let Some(result) = rhtslib_reimplements::rhtslib_read(htsfile, itr, &mut record) {
        // Should we return an Error here or just continue to the next iteration?
        result.with_context(|| {
            format!(
                "Encountered faulty read for {}",
                tr_region.reference_info.get_fetch_definition_s()
            )
        })?;

        if record.is_duplicate() || record.is_supplementary() || record.is_quality_check_failed() {
            // Ignore duplicate, supplementary, low quality reads
            continue;
        }
        if record.pos() >= tr_region.reference_info.start - flank as i64
            || record.reference_end() <= tr_region.reference_info.end + flank as i64
        {
            // This is not an enclosing read: skip
            continue;
        }

        let starting_pos = record.pos(); // 0-based
        let tr_region_len = allele_length_from_cigar(
            &record.cigar(),
            starting_pos,
            tr_region.reference_info.start,
            tr_region.reference_info.end,
        )?;
        if tr_region_len % tr_region.reference_info.period != 0 {
            // TR length is not a multiple of period: skip
            continue;
        }
        let tr_len = tr_region_len / tr_region.reference_info.period;

        allele_lengths
            .entry(tr_len)
            .and_modify(|counter| *counter += 1.0)
            .or_insert(1.0);
    }
    if !allele_lengths.is_empty() {
        tr_region.set_allele_lengths(Some(allele_lengths));
    }

    Ok(())
}

fn allele_length_from_cigar(
    cigar: &CigarStringView,
    mut current_pos: i64,
    tr_start: i64,
    tr_end: i64,
) -> Result<i64> {
    let mut tr_region_len = 0;

    for op in cigar {
        let consumes_r = cigar::consumes_ref(op);
        let advances_tr = cigar::advances_tr_len(op);
        let len = i64::from(op.len());

        if consumes_r && !advances_tr {
            current_pos += len;
        } else if !consumes_r && advances_tr {
            // We can be sure that current_pos <= tr_region.reference_info.end since this is tested at
            // the end of every loop, we don't need to test here again
            if current_pos >= tr_start {
                tr_region_len += len;
            }
        } else if consumes_r && advances_tr {
            // BAM and BED coordinate systems are half-open, utils::range_overlap assumes closed => decrement end positions
            tr_region_len +=
                utils::range_overlap(current_pos, current_pos + len - 1, tr_start, tr_end - 1)?;
            current_pos += len;
        }

        if current_pos >= tr_end {
            break;
        }
    }

    Ok(tr_region_len)
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::record::{Cigar, CigarString};

    use super::*;

    #[test]
    fn tr_length_from_cigar_match() {
        let cigar_start = 20;
        let cigar = CigarString(vec![Cigar::Match(100)]).into_view(cigar_start);
        let tr_region_len = allele_length_from_cigar(&cigar, cigar_start, 40, 50).unwrap();

        assert_eq!(10, tr_region_len);
    }
    #[test]
    fn tr_length_from_cigar_ins() {
        let cigar_start = 20;
        let cigar = CigarString(vec![Cigar::Match(20), Cigar::Ins(6), Cigar::Match(54)])
            .into_view(cigar_start);
        let tr_region_len = allele_length_from_cigar(&cigar, cigar_start, 40, 50).unwrap();

        assert_eq!(16, tr_region_len);
    }
    #[test]
    fn tr_length_from_cigar_del() {
        let cigar_start = 20;
        let cigar = CigarString(vec![Cigar::Match(20), Cigar::Del(5), Cigar::Match(54)])
            .into_view(cigar_start);
        let tr_region_len = allele_length_from_cigar(&cigar, cigar_start, 40, 50).unwrap();

        assert_eq!(5, tr_region_len);
    }
}
