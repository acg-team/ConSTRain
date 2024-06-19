//! # Utils
//!
//! Modules containing utility functions and structs for the ConSTRain library.
//! This top-level module contains miscellaneous utility functions,
//! the sub-modules contain functions related to specific functionality.
pub mod cigar_utils;
use ndarray::prelude::*;
pub mod io_utils;
use std::{cmp, error::Error, path::Path};

/// A structure to represent Copy Number Variants
#[derive(Debug, serde::Deserialize)]
pub struct CopyNumberVariant {
    pub seqname: String,
    // CNV struct follows 0-based half-open coordinate system: [start, end)
    pub start: i64,
    pub end: i64,
    pub cn: usize,
}

/// Determine the overlap between two ranges, each specified by their start
/// and end coordinates.
///
/// # Examples
///
/// ```
/// let a: Vec<i64> = vec![10, 15];
/// let b: Vec<i64> = vec![13, 25];
/// let overlap = constrain::utils::range_overlap(a[0], a[1], b[0], b[1]);
/// assert_eq!(3, overlap);
/// ```
pub fn range_overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> i64 {
    cmp::max(0, cmp::min(a_end, b_end) - cmp::max(a_start, b_start) + 1)
}

/// Infer a sample name from the filepath of an alignment file
///
/// # Examples
///
/// ```
/// let filepath = "./path/to/alignment.bam";
/// let sample_name = constrain::utils::sample_name_from_path(filepath).unwrap();
///
/// assert_eq!("alignment", sample_name);
/// ```
pub fn sample_name_from_path(filepath: &str) -> Result<String, Box<dyn Error>> {
    let name = Path::new(filepath)
        .file_stem()
        .ok_or("Could not infer sample name from path")?
        .to_str()
        .ok_or("Could not infer sample name from path")?;

    Ok(String::from(name))
}

/// Zero pad array `a` to the given `min_len`.
/// If `a` is already of `min_len` or longer, do nothing and return  `a`
///
/// # Examples
/// ```
/// use constrain::utils::zero_pad_if_shorter;
/// use ndarray::prelude::*;
///
/// let left = arr1(&[13., 12., 11.]);
/// let right = arr1(&[13., 12., 11.]);
/// assert_eq!(left, zero_pad_if_shorter(right, 3));
///
/// let left = arr1(&[13., 12., 11., 0., 0.]);
/// let right = arr1(&[13., 12., 11.]);
/// assert_eq!(left, zero_pad_if_shorter(right, 5));
/// ```
pub fn zero_pad_if_shorter(
    a: Array<f32, Dim<[usize; 1]>>,
    min_len: usize,
) -> Array<f32, Dim<[usize; 1]>> {
    if a.len() >= min_len {
        return a;
    }
    let mut padded_a = Array::<f32, _>::zeros(min_len);
    padded_a.slice_mut(s![..a.len()]).assign(&a);
    padded_a
}

/// Number of partitions that exist for integers 0 - 50 (see [https://oeis.org/A000041](https://oeis.org/A000041)).
/// For 0 - 1000, see here: [https://oeis.org/A000041/b000041.txt](https://oeis.org/A000041/b000041.txt)
pub const N_PARTITIONS: &[usize] = &[
    1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792,
    1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977,
    21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525,
    204226,
];

mod tests {
    use super::*;

    #[test]
    fn zero_padding() {
        let left = arr1(&[13., 12., 11.]);
        let right = arr1(&[13., 12., 11.]);
        assert_eq!(left, zero_pad_if_shorter(right, 3));

        let left = arr1(&[13., 12., 11., 0., 0.]);
        let right = arr1(&[13., 12., 11.]);
        assert_eq!(left, zero_pad_if_shorter(right, 5));
    }
}
