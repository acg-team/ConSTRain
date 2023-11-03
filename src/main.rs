use std::{cmp, fs, io, collections::HashMap, env};

use rayon::prelude::*;
use rust_htslib::{bam, bam::{Read, ext::BamRecordExtensions}, bam::record::Cigar};
use noodles::bed;

fn main() {
    // USAGE: ./target/debug/CN-guided-STR-genotying data/repeats/APC_repeats.tsv data/alignments/patient_3.bam 4
    // USAGE: CN-guided-STR-genotying <bed file> <bam file> <threads>
    // Eventually make a nicer CLI using Clap crate?
    let args: Vec<String> = env::args().collect();
    let bed_path = &args[1];
    let bam_path = &args[2];
    let n_threads = &args[3].parse::<usize>().unwrap();

    rayon::ThreadPoolBuilder::new().num_threads(*n_threads).build_global().unwrap();
    let mut reader = fs::File::open(bed_path)
        .map(io::BufReader::new)
        .map(bed::Reader::new)
        .unwrap();

    let mut str_regions = Vec::new();
    for result in reader.records::<3>() {
        str_regions.push(result.unwrap());
    }

    let flank = 4;
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
                let consumes_r = cigar_consumes_ref(op);
                let advances_str = cigar_advances_str_len(op);
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
                    str_len += range_overlap(current_pos, current_pos + len - 1, str_region.1, str_region.2);
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

fn cigar_consumes_ref(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Del(_)       |
        Cigar::RefSkip(_)   |
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

fn cigar_consumes_query(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Ins(_)       |
        Cigar::SoftClip(_)  | 
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

fn cigar_advances_str_len(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Ins(_)       |
        // Cigar::SoftClip(_)  | // softclip consumes reference but shouldn't advance STR length (right?)
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

fn range_overlap(left_start: i64, left_end: i64, right_start: i64, right_end: i64) -> i64 {
    cmp::max(0, cmp::min(left_end, right_end) - cmp::max(left_start, right_start) + 1 )
}
