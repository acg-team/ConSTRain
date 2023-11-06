use std::cmp;

use rust_htslib::bam::record::Cigar;

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

pub fn cigar_advances_str_len(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Ins(_)       |
        // Cigar::SoftClip(_)  | // softclip consumes reference but shouldn't advance STR length (right?)
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}

pub fn range_overlap(left_start: i64, left_end: i64, right_start: i64, right_end: i64) -> i64 {
    cmp::max(0, cmp::min(left_end, right_end) - cmp::max(left_start, right_start) + 1 )
}