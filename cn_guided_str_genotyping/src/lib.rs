pub mod utils;

use noodles::bed;
use rayon::prelude::*;
use rust_htslib::bam::{self, {Read, ext::BamRecordExtensions}};
use std::{fs::File, io::BufReader, collections::HashMap};

pub fn alelle_lengths_from_cigar(bed_path: &String, bam_path: &String, flank: i64) {
    let mut reader = File::open(bed_path)
        .map(BufReader::new)
        .map(bed::Reader::new)
        .unwrap();
    
    let mut str_regions = Vec::new();
    for result in reader.records::<3>() {
        str_regions.push(result.unwrap());
    }
    
    let mut results = Vec::new();
    str_regions.par_iter().map(|str_region| {
        let str_region: (&str, i64, i64) = (
            str_region.reference_sequence_name(), 
            str_region.start_position().get() as i64 - 2, 
            str_region.end_position().get() as i64 - 1
        );
        let mut alelle_lengths: HashMap<i64, usize> = HashMap::new();

        let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
        bam.fetch(str_region).unwrap();

        for record in bam.records() {
            let record = record.unwrap();
            if record.is_duplicate() || record.is_supplementary() {
                // Ignore duplicate or supplementary reads
                continue
            }

            if record.pos() >= str_region.1 - flank || record.reference_end() <= str_region.2 + flank {
                // This is not an enclosing read: skip
                continue
            }
            
            let mut current_pos = record.pos();
            let mut str_len = 0;
            for op in &record.cigar() {
                let consumes_r = utils::cigar_consumes_ref(op);
                let advances_str = utils::cigar_advances_str_len(op);
                let len = op.len() as i64;

                if consumes_r && !advances_str {
                    current_pos += len;                 
                } else if !consumes_r & advances_str {
                    // We can be sure that current_pos <= str_region.2 since this is tested at
                    // the end of every loop, we don't need to test here again
                    if current_pos >= str_region.1 {
                        str_len += len;
                    }
                } else if consumes_r && advances_str {                
                    str_len += utils::range_overlap(current_pos, current_pos + len - 1, str_region.1, str_region.2);
                    current_pos += len;
                }

                if current_pos >= str_region.2 {
                    break
                }
            }
            alelle_lengths
                .entry(str_len)
                .and_modify(|counter| *counter += 1)
                .or_insert(1);            
        }
        (str_region, alelle_lengths)
    }).collect_into_vec(&mut results);

    for r in results {
        println!("{:?}", r);
    }
}
