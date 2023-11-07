mod utils;

use crate::utils::{ReferenceRepeat, io_utils};

use rayon::prelude::*;
use rust_htslib::bam::{self, {Read, ext::BamRecordExtensions}};
use std::collections::HashMap;

pub fn call_tr_allele_lengths(bed_path: &String, bam_path: &String, flank: i64) {
    let mut tr_regions = Vec::new();
    io_utils::trs_from_bed(bed_path, &mut tr_regions);
    
    let mut results: Vec<(&ReferenceRepeat, HashMap<i64, usize>)> = Vec::new();
    tr_regions.par_iter().map(|tr_region| {
        alelle_lengths_from_cigar(tr_region, bam_path, flank)
    }).collect_into_vec(&mut results);

    for r in results {
        println!("{:?}", r);
    }
}

fn alelle_lengths_from_cigar<'a> (
    tr_region  : &'a ReferenceRepeat, 
    bam_path    : &String, 
    flank       : i64
) -> (&'a ReferenceRepeat, HashMap<i64, usize>) {
    let fetch_req = tr_region.get_fetch_definition();
    let mut alelle_lengths: HashMap<i64, usize> = HashMap::new();

    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    bam.fetch(fetch_req).unwrap();

    for record in bam.records() {
        let record = record.unwrap();
        if record.is_duplicate() || record.is_supplementary() {
            // Ignore duplicate or supplementary reads
            continue
        }

        if record.pos() >= tr_region.start - flank || record.reference_end() <= tr_region.end + flank {
            // This is not an enclosing read: skip
            continue
        }
        
        let mut current_pos = record.pos();
        let mut tr_len = 0;
        for op in &record.cigar() {
            let consumes_r = utils::cigar_consumes_ref(op);
            let advances_tr = utils::cigar_advances_tr_len(op);
            let len = op.len() as i64;

            if consumes_r && !advances_tr {
                current_pos += len;                 
            } else if !consumes_r & advances_tr {
                // We can be sure that current_pos <= tr_region.end since this is tested at
                // the end of every loop, we don't need to test here again
                if current_pos >= tr_region.start {
                    tr_len += len;
                }
            } else if consumes_r && advances_tr {        
                // BAM and BED coordinate system is half-open, utils::range_overlap assumes closed => decrement end positions         
                tr_len += utils::range_overlap(current_pos, current_pos + len - 1, tr_region.start, tr_region.end - 1);

                current_pos += len;
            }

            if current_pos >= tr_region.end {
                break
            }
        }
        if tr_len % tr_region.period != 0 {
            // TR length is not a multiple of period: skip
            continue;
        }
        tr_len /= tr_region.period;

        alelle_lengths
            .entry(tr_len)
            .and_modify(|counter| *counter += 1)
            .or_insert(1);            
    }
    (tr_region, alelle_lengths)
}
