pub mod genotyping;
pub mod repeat;
pub mod utils;

use crate::repeat::TandemRepeat;
use crate::utils::{cigar_utils, CopyNumberVariant};

use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};
use std::collections::HashMap;

pub fn alelle_lengths_from_cigar(tr_region: &mut TandemRepeat, bam_path: &str, flank: i64) {
    let fetch_req = tr_region.reference_info.get_fetch_definition();
    let mut allele_lengths: HashMap<i64, f32> = HashMap::new();

    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    bam.fetch(fetch_req).unwrap();

    for record in bam.records() {
        let record = record.unwrap();

        if record.is_duplicate() || record.is_supplementary() || record.is_quality_check_failed() {
            // Ignore duplicate, supplementary, low quality reads
            continue;
        }
        if record.pos() >= tr_region.reference_info.start - flank
            || record.reference_end() <= tr_region.reference_info.end + flank
        {
            // This is not an enclosing read: skip
            continue;
        }

        let mut current_pos = record.pos();
        let mut tr_len = 0;
        for op in &record.cigar() {
            let consumes_r = cigar_utils::cigar_consumes_ref(op);
            let advances_tr = cigar_utils::cigar_advances_tr_len(op);
            let len = op.len() as i64;

            if consumes_r && !advances_tr {
                current_pos += len;
            } else if !consumes_r && advances_tr {
                // We can be sure that current_pos <= tr_region.reference_info.end since this is tested at
                // the end of every loop, we don't need to test here again
                if current_pos >= tr_region.reference_info.start {
                    tr_len += len;
                }
            } else if consumes_r && advances_tr {
                // BAM and BED coordinate systems are half-open, utils::range_overlap assumes closed => decrement end positions
                tr_len += utils::range_overlap(
                    current_pos,
                    current_pos + len - 1,
                    tr_region.reference_info.start,
                    tr_region.reference_info.end - 1,
                );
                current_pos += len;
            }

            if current_pos >= tr_region.reference_info.end {
                break;
            }
        }
        if tr_len % tr_region.reference_info.period != 0 {
            // TR length is not a multiple of period: skip
            continue;
        }
        tr_len /= tr_region.reference_info.period;

        // TODO: Should try to do this directly on the TandemRepeat struct
        allele_lengths
            .entry(tr_len)
            .and_modify(|counter| *counter += 1.0)
            .or_insert(1.0);
        tr_region.n_mapped_reads += 1;
    }

    if !allele_lengths.is_empty() {
        tr_region.set_allele_lengths(Some(allele_lengths));
    }
}


pub fn tr_cn_from_cnvs(tr_region: &mut TandemRepeat, cnv_regions: &[CopyNumberVariant]) {
    // Currently extremely basic one vs all comparison strategy
    // Check how GenomicRanges R library (or bedtools?) finds range
    // overlaps between two lists of entities
    for cnv in cnv_regions {
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
            // TR partially overlaps CNV, impossible to set sensible CN for TR, so set to 0 so it gets skipped
            // Not strictly correct, maybe add this concept explicitly to TandemRepeat
            tr_region.set_cn(0);
            return;
        }
    }
}
