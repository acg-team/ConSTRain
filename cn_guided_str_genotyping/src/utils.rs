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

// Number of partitions that exist for integers 0 - 50 (see https://oeis.org/A000041)
pub const N_PARTITIONS: &'static [usize] = &[
    1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792,
    1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977,
    21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525,
    204226,
];
