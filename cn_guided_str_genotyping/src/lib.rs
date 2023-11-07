mod utils;
mod repeat;

use crate::repeat::TandemRepeat;
use crate::utils::io_utils;

use rayon::prelude::*;
use rust_htslib::bam::{self, {Read, ext::BamRecordExtensions}};
use std::collections::HashMap;

pub fn call_tr_allele_lengths(bed_path: &String, bam_path: &String, ploidy_path: &String, flank: i64) {
    let mut tr_regions = Vec::new();
    io_utils::trs_from_bed(bed_path, ploidy_path, &mut tr_regions);

    // Here, check if segment CNV data is supplied. If so: apply CNV calls to STR
    // copy_number values. Need to check for overlap between STRs and CNVs => check GenomicRanges R library
    
    let mut results: Vec<&TandemRepeat> = Vec::new();
    tr_regions.par_iter_mut().map(|tr_region| {
        if tr_region.copy_number != 0 {
            alelle_lengths_from_cigar(tr_region, bam_path, flank)
        } else {
            tr_region
        }
    }).collect_into_vec(&mut results);

    for r in results {
        println!("{:?}", r);
    }
}

fn alelle_lengths_from_cigar<'a> (
    tr_region   : &'a mut TandemRepeat, 
    bam_path    : &String, 
    flank       : i64
) -> &'a TandemRepeat {
    let fetch_req = tr_region.reference_info.get_fetch_definition();
    let mut alelle_lengths: HashMap<i64, usize> = HashMap::new();

    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    bam.fetch(fetch_req).unwrap();

    for record in bam.records() {
        let record = record.unwrap();
        if record.is_duplicate() || record.is_supplementary() {
            // Ignore duplicate or supplementary reads
            continue
        }

        if record.pos() >= tr_region.reference_info.start - flank || record.reference_end() <= tr_region.reference_info.end + flank {
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
                // We can be sure that current_pos <= tr_region.reference_info.end since this is tested at
                // the end of every loop, we don't need to test here again
                if current_pos >= tr_region.reference_info.start {
                    tr_len += len;
                }
            } else if consumes_r && advances_tr {        
                // BAM and BED coordinate system is half-open, utils::range_overlap assumes closed => decrement end positions         
                tr_len += utils::range_overlap(current_pos, current_pos + len - 1, tr_region.reference_info.start, tr_region.reference_info.end - 1);

                current_pos += len;
            }

            if current_pos >= tr_region.reference_info.end {
                break
            }
        }
        if tr_len % tr_region.reference_info.period != 0 {
            // TR length is not a multiple of period: skip
            continue;
        }
        tr_len /= tr_region.reference_info.period;

        alelle_lengths
            .entry(tr_len)
            .and_modify(|counter| *counter += 1)
            .or_insert(1);            
    }

    tr_region.set_genotyped();
    if !alelle_lengths.is_empty() {        
        tr_region.set_allele_lengths(Some(alelle_lengths));
    }    

    tr_region
}
