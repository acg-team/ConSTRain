pub mod genotyping;
pub mod repeat;
pub mod utils;

use crate::repeat::TandemRepeat;
use crate::utils::{cigar_utils, CopyNumberVariant};

use genotyping::partitions;
use ndarray::{Array, Dim};
use rust_htslib::{
    bam::{ext::BamRecordExtensions, record::CigarStringView, Record},
    errors::Error as htsError,
    htslib,
};
use std::collections::HashMap;

pub fn fetch_allele_lengths(
    tr_region: &mut TandemRepeat,
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    flank: usize,
) -> Result<(), String> {
    let mut allele_lengths: HashMap<i64, f32> = HashMap::new();
    let mut record = Record::new();
    while let Some(result) = rhtslib_read(htsfile, itr, &mut record) {
        if result.is_err() {            
            eprintln!(
                "Faulty read for region {}",
                tr_region.reference_info.get_fetch_definition_s()
            );
            continue;
        }

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

        let starting_pos = record.pos();
        let tr_region_len = allele_length_from_cigar(
            &record.cigar(),
            starting_pos,
            tr_region.reference_info.start,
            tr_region.reference_info.end,
        );
        if tr_region_len % tr_region.reference_info.period != 0 {
            // TR length is not a multiple of period: skip
            continue;
        }
        let tr_len = tr_region_len / tr_region.reference_info.period;
        // TODO: Should try to do this directly on the TandemRepeat struct
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
) -> i64 {
    let mut tr_region_len = 0;
    for op in cigar {
        let consumes_r = cigar_utils::cigar_consumes_ref(op);
        let advances_tr = cigar_utils::cigar_advances_tr_len(op);
        let len = op.len() as i64;

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
                utils::range_overlap(current_pos, current_pos + len - 1, tr_start, tr_end - 1);
            current_pos += len;
        }

        if current_pos >= tr_end {
            break;
        }
    }

    tr_region_len
}

pub fn tr_cn_from_cnvs(tr_region: &mut TandemRepeat, cnv_regions: &[CopyNumberVariant]) {
    // Currently extremely basic one vs all comparison strategy
    // Check how GenomicRanges R library (or bedtools?) finds range
    // overlaps between two lists of entities
    for cnv in cnv_regions {
        if tr_region.reference_info.seqname != cnv.seqname {
            continue;
        }
        let overlap = utils::range_overlap(
            tr_region.reference_info.start,
            tr_region.reference_info.end - 1,
            cnv.start,
            cnv.end - 1,
        );
        if overlap == tr_region.reference_info.end - tr_region.reference_info.start {
            tr_region.set_cn(cnv.cn);
            return;
        } else if overlap > 0 {
            // TR partially overlaps CNV, impossible to set sensible CN for TR. Set to 0 so it gets skipped
            // Not strictly correct, maybe add this concept explicitly to TandemRepeat
            tr_region.set_cn(0);
            return;
        }
    }
}

pub fn make_partitions_map(copy_numbers: &[usize]) -> HashMap<usize, Array<f32, Dim<[usize; 2]>>> {
    let mut map: HashMap<usize, Array<f32, Dim<[usize; 2]>>> = HashMap::new();
    for cn in copy_numbers {
        if *cn == 0 {
            continue;
        }
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
        let tr_region_len = allele_length_from_cigar(&cigar, cigar_start, 40, 50);

        assert_eq!(10, tr_region_len);
    }
    #[test]
    fn tr_length_from_cigar_ins() {
        let cigar_start = 20;
        let cigar = CigarString(vec![Cigar::Match(20), Cigar::Ins(6), Cigar::Match(54)])
            .into_view(cigar_start);
        let tr_region_len = allele_length_from_cigar(&cigar, cigar_start, 40, 50);

        assert_eq!(16, tr_region_len);
    }
    #[test]
    fn tr_length_from_cigar_del() {
        let cigar_start = 20;
        let cigar = CigarString(vec![Cigar::Match(20), Cigar::Del(5), Cigar::Match(54)])
            .into_view(cigar_start);
        let tr_region_len = allele_length_from_cigar(&cigar, cigar_start, 40, 50);

        assert_eq!(5, tr_region_len);
    }
}

// Functions below that are prefixed with 'rhtslib_' are private functions in
// rust_htslib, and are copied from there to be used in this library.
// It would be nicer to use rust_htslib's structs (e.g., IndexedReader) for reading
// alignment files, but those structs caused segfaults when reading CRAM files (not BAM, interestingly)
// when they were dropped, even on trivial tests. -> submit issue on GitHub?
fn rhtslib_read(
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    record: &mut Record,
) -> Option<Result<(), htsError>> {
    match rhtslib_itr_next(htsfile, itr, &mut record.inner as *mut htslib::bam1_t) {
        -1 => None,
        -2 => Some(Err(htsError::BamTruncatedRecord)),
        -4 => Some(Err(htsError::BamInvalidRecord)),
        _ => {
            // record.set_header(Rc::clone(&self.header)); // this does not seem to be necessary?
            Some(Ok(()))
        }
    }
}

fn rhtslib_itr_next(
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    record: *mut htslib::bam1_t,
) -> i32 {
    unsafe {
        htslib::hts_itr_next(
            (*htsfile).fp.bgzf,
            itr,
            record as *mut ::std::os::raw::c_void,
            htsfile as *mut ::std::os::raw::c_void,
        )
    }
}
