pub mod genotyping;
pub mod repeat;
pub mod rhtslib_reimplement;
pub mod utils;

use crate::{
    repeat::TandemRepeat,
    utils::{cigar_utils, CopyNumberVariant},
};

use anyhow::{Context, Result};
use genotyping::partitions;
use ndarray::{Array, Dim};
use rust_htslib::{
    bam::{ext::BamRecordExtensions, record::CigarStringView, Record},
    htslib::{self, htsFile},
};
use std::{collections::HashMap, ffi, sync::Arc};

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
    cnv_regions: &[CopyNumberVariant],
    partitions_map: &Arc<HashMap<usize, Array<f32, Dim<[usize; 2]>>>>,
    alignment: &str,
    reference: Option<&String>,
    flanksize: usize,
    reads_per_allele: usize,
    tidx: usize,
) {
    eprintln!("Launching thread {tidx}");

    let (htsfile, idx, header) =
        thread_setup(alignment, reference).expect("Error during thread setup");

    for tr_region in tr_regions {
        let fetch_request = tr_region.reference_info.get_fetch_definition_s();
        if !cnv_regions.is_empty() {
            if let Err(e) = tr_cn_from_cnvs(tr_region, &cnv_regions) {
                eprintln!("Thread {tidx}: Error setting copy number, skipping locus {fetch_request}: {e:?}");
                continue;
            };
        }

        let Ok(itr) =
            rhtslib_reimplement::rhtslib_fetch_by_str(idx, header, fetch_request.as_bytes())
        else {
            eprintln!("Thread {tidx}: Error fetching reads, skipping locus {fetch_request}");
            continue;
        };

        if let Err(e) = extract_allele_lengths(tr_region, htsfile, itr, flanksize) {
            eprintln!("Thread {tidx}: Error extracting allele lengths, skipping locus {fetch_request}: {e:?}");
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

        if let Err(e) =
            // Should we do Arc::clone for each tr_region or could we do it once per thread instead?
            genotyping::estimate_genotype(
                tr_region,
                reads_per_allele,
                Arc::clone(&partitions_map),
            )
        {
            eprintln!(
                "Thread {tidx}: Could not estimate genotype for locus {fetch_request}: {e:?}"
            );
            continue;
        }
    }
    unsafe {
        htslib::hts_close(htsfile);
    }

    eprintln!("Finished on thread {tidx}");
}

fn thread_setup(
    alignment_path: &str,
    reference: Option<&String>,
) -> Result<(*mut htsFile, *mut htslib::hts_idx_t, *mut htslib::sam_hdr_t)> {
    let htsfile = rhtslib_reimplement::rhtslib_from_path(alignment_path)?;
    let header: *mut htslib::sam_hdr_t = unsafe { htslib::sam_hdr_read(htsfile) };
    let c_str = ffi::CString::new(alignment_path)
        .context("Internal 0 byte contained in alignment file name")?;
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
            rhtslib_reimplement::rhtslib_set_reference(htsfile, reference)?;
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
    while let Some(result) = rhtslib_reimplement::rhtslib_read(htsfile, itr, &mut record) {
        // Should we return an Error here or just continue to the next iteration?
        result.with_context(|| {
            format!(
                "Encountered faulty read for {}",
                tr_region.reference_info.get_fetch_definition_s()
            )
        })?;
        // if result.is_err() {
        //     // eprintln!(
        //     //     "Faulty read for region {}",
        //     //     tr_region.reference_info.get_fetch_definition_s()
        //     // );
        //     // continue;
        // }

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
        // TODO: Could try to do this directly on the TandemRepeat struct
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
        let consumes_r = cigar_utils::cigar_consumes_ref(op);
        let advances_tr = cigar_utils::cigar_advances_tr_len(op);
        // let len = op.len() as i64;
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

pub fn tr_cn_from_cnvs(
    tr_region: &mut TandemRepeat,
    cnv_regions: &[CopyNumberVariant],
) -> Result<()> {
    // Currently extremely basic one vs all comparison strategy
    // Check how GenomicRanges R library (or bedtools?) finds range
    // overlaps between two lists of entities
    let region_len = tr_region.reference_info.end - tr_region.reference_info.start;
    for cnv in cnv_regions {
        if tr_region.reference_info.seqname != cnv.seqname {
            continue;
        }
        let overlap = utils::range_overlap(
            tr_region.reference_info.start,
            tr_region.reference_info.end - 1,
            cnv.start,
            cnv.end - 1,
        )?;

        if overlap == region_len {
            tr_region.set_cn(cnv.cn);
            // TRs can intersect with at most one CNV, we found a hit so we can return
            return Ok(());
        } else if overlap > 0 {
            // TR partially overlaps CNV, impossible to set sensible CN for TR. Set to 0 so it gets skipped
            // Should we return an Err here instead?
            tr_region.set_cn(0);
            // TRs can intersect with at most one CNV, we found a hit so we can return
            return Ok(());
        }
    }
    Ok(())
}

pub fn make_partitions_map(copy_numbers: &[usize]) -> HashMap<usize, Array<f32, Dim<[usize; 2]>>> {
    let mut map: HashMap<usize, Array<f32, Dim<[usize; 2]>>> = HashMap::new();
    for cn in copy_numbers {
        // if *cn == 0 {
        //     continue;
        // }
        map.insert(*cn, partitions(*cn));
    }

    map
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
