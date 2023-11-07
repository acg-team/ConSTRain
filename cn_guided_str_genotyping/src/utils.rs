pub mod io_utils;

use std::{cmp, collections::HashMap};

use rust_htslib::bam::record::Cigar;

#[derive(Debug, serde::Deserialize)]
pub struct ReferenceRepeat {
    pub seqname: String,
    // Tandem repeat struct follows 0-based half-open coordinate system: [start, end)
    pub start: i64,
    pub end: i64,
    pub period: i64,
    pub unit: String,
}

impl ReferenceRepeat {
    pub fn get_fetch_definition(&self) -> (&str, i64, i64) {
        (self.seqname.as_str(), self.start, self.end)
    }
}

#[derive(Debug)]
pub struct TandemRepeat<'a> {
    pub reference_info: &'a ReferenceRepeat,
    pub is_genotyped: bool,
    pub allele_lengths: Option<HashMap<i64, usize>>,
}

pub fn cigar_consumes_ref(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Del(_)       |
        Cigar::RefSkip(_)   |
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

pub fn cigar_consumes_query(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Ins(_)       |
        Cigar::SoftClip(_)  | 
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

pub fn cigar_advances_tr_len(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Ins(_)       |
        // Cigar::SoftClip(_)  | // softclip consumes reference but shouldn't advance TR length (right?)
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

pub fn range_overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> i64 {
    cmp::max(0, cmp::min(a_end, b_end) - cmp::max(a_start, b_start) + 1 )
}