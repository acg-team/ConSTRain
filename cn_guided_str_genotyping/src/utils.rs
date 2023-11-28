pub mod cigar_utils;
pub mod io_utils;

use std::cmp;

#[derive(Debug, serde::Deserialize)]
pub struct CopyNumberVariant {
    pub seqname: String,
    // CNV struct follows 0-based half-open coordinate system: [start, end)
    pub start: i64,
    pub end: i64,
    pub cn: usize,
}

pub fn range_overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> i64 {
    cmp::max(0, cmp::min(a_end, b_end) - cmp::max(a_start, b_start) + 1)
}
