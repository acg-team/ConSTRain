//! # Estimating genotypes from allele length distributions
use anyhow::{bail, Context, Result};
use ndarray::prelude::*;
use std::{cmp::Ordering, collections::HashMap, sync::Arc};

use crate::{
    repeat::TandemRepeat,
    utils::{self, N_PARTITIONS},
};

/// Estimate the most likely underlying genotype that produced the
/// observed allele distribution stored on the `tr_region` struct.
pub fn estimate_genotype(
    tr_region: &mut TandemRepeat,
    min_reads_per_allele: usize,
    partitions_map: Arc<HashMap<usize, Array<f32, Dim<[usize; 2]>>>>,
) -> Result<()> {
    // Copy number for locus is 0, there is no need to estimate
    if tr_region.copy_number == 0 {
        return Ok(());
    }
    // If no reads were mapped to this locus, there is no need to estimate
    let Some(n_mapped_reads) = tr_region.get_n_mapped_reads() else {
        return Ok(());
    };
    if n_mapped_reads < min_reads_per_allele * tr_region.copy_number {
        // Not enough reads were mapped to this locus, refuse to estimate
        return Ok(());
    }

    let (allele_lengths, mut counts) = tr_region.allele_freqs_as_ndarrays(Some("freq"));
    // might need this error_constant at some point if trying to infer CN from read distribution
    // let mut error_constant = 0.;
    match counts.len().cmp(&tr_region.copy_number) {
        Ordering::Less => counts = utils::zero_pad_if_shorter(counts, tr_region.copy_number),
        Ordering::Equal => (),
        Ordering::Greater => counts = counts.slice_move(s![..tr_region.copy_number]), // and maybe set error_constant?
    };

    let partitions = partitions_map
        .get(&tr_region.copy_number)
        .with_context(|| {
            format!(
                "Repeat {} has copy number {}, for which no partitions exist",
                tr_region.reference_info.get_fetch_definition_s(),
                tr_region.copy_number
            )
        })?;

    // At least `threshold_val` number of reads must support a give allele length for it to be considered.
    // If the observed count for a specific allele is closer to 0 than to the expected count for an allele
    // that is present once, we ignore this allele.
    let threshold_val = n_mapped_reads as f32 / tr_region.copy_number as f32 * 0.5;
    let valid_partition_idxs = find_valid_partition_idxs(partitions, &counts, threshold_val);

    if valid_partition_idxs.is_empty() {
        // Not a single allele length observed with frequency over threshold, refuse to estimate
        return Ok(());
    } else if valid_partition_idxs.len() == 1 {
        // Only one allele length observed with frequency over threshold, we can exit early
        let genotype: Vec<(i64, f32)> = vec![(allele_lengths[0], tr_region.copy_number as f32)];
        tr_region.genotype = Some(genotype);
        return Ok(());
    }

    let valid_partitions = partitions.select(Axis(0), &valid_partition_idxs);

    let argmin = most_likely_partition_idx(
        &valid_partitions,
        n_mapped_reads,
        tr_region.copy_number,
        &counts,
    )?;
    let most_likely_partition = valid_partitions.slice(s![argmin, ..]);

    let mut genotype: Vec<(i64, f32)> = allele_lengths
        .iter()
        .zip(most_likely_partition.iter())
        .filter_map(|(x, y)| if *y > 0. { Some((*x, *y)) } else { None })
        .collect();
    genotype.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    tr_region.genotype = Some(genotype);

    Ok(())
}

/// Not all partitions may be relevant for estimating the genotype of the
/// observed allele length distribution. Specifically, some partitions may contain
/// more distinct alleles than have been observed in the alignment. For example, we
/// don't want to consider partition `[1., 1., 1.]` for a TR with copy number 3 when
/// the observed allele distribution is `[20., 8., 0.]`
fn find_valid_partition_idxs(
    partitions: &Array<f32, Dim<[usize; 2]>>,
    counts: &Array<f32, Dim<[usize; 1]>>,
    threshold_value: f32,
) -> Vec<usize> {
    // First, find how many alleles are observed that are supported by
    // enough reads (i.e., have more reads than `threshold_value`)
    let mut threshold_idx = 0;
    for count in counts {
        if *count > threshold_value {
            threshold_idx += 1;
        } else {
            break;
        }
    }

    // Next, find partitions where the first `threshold_idx` positions sum to
    // the number of alleles needed to make the current copy number
    let valid_partitions = partitions
        .slice(s![.., ..threshold_idx])
        .sum_axis(Axis(1))
        .iter()
        .enumerate()
        .filter_map(|(idx, val)| {
            if *val as usize == partitions.shape()[1] {
                Some(idx)
            } else {
                None
            }
        })
        .collect();

    valid_partitions
}

/// Return the indexes of the most likely genotype to underly the observed allele distribution
fn most_likely_partition_idx(
    partitions: &Array<f32, Dim<[usize; 2]>>,
    n_mapped_reads: usize,
    copy_number: usize,
    counts: &Array<f32, Dim<[usize; 1]>>,
) -> Result<usize> {
    let mut errors = partitions.mapv(|a| a * (n_mapped_reads / copy_number) as f32);
    errors = errors - counts.slice(s![..copy_number]);
    errors.mapv_inplace(|x| x.powi(2));

    let error_sums = errors.sum_axis(Axis(1));
    // we can unwrap here: `min_by()` returns None if iterator is empty and we check earlier that it's not
    let min = *error_sums.iter().min_by(|a, b| a.total_cmp(b)).unwrap();

    let argmin: Vec<usize> = error_sums
        .iter()
        .enumerate()
        .filter_map(|(x, y)| {
            // For each value, check if it equals the minimum
            if (*y - min).abs() < f32::EPSILON {
                Some(x)
            } else {
                None
            }
        })
        .collect();

    match argmin.len().cmp(&1) {
        Ordering::Less => {
            bail!("No most likely allele distribution was found. This should not happen")
        }
        Ordering::Equal => Ok(argmin[0]),
        Ordering::Greater => bail!("Multiple allele distributions are equally likely"),
    }
}

/// Generate all partitions of for integer 'n' ([Wikipedia](https://en.wikipedia.org/wiki/Integer_partition), [Mathworld](https://mathworld.wolfram.com/Partition.html)).
/// All partitions are padded with zeroes to be of size 'n', and values are represented by f32s.
/// This is because the partitions will be used to calculate a mean-squared errors later.
/// Implementation based on the iterative algorithm described in <https://jeromekelleher.net/category/combinatorics>.
pub fn partitions(n: usize) -> Array<f32, Dim<[usize; 2]>> {
    if n == 0 {
        return arr2(&[[]]);
    }
    let mut results = Array::zeros((N_PARTITIONS[n], n));
    let mut results_idx = N_PARTITIONS[n] - 1;

    let mut a = vec![0; n + 1];
    let mut k = 1;
    a[1] = n;
    while k != 0 {
        let x = a[k - 1] + 1;
        let mut y = a[k] - 1;
        k -= 1;
        while x <= y {
            a[k] = x;
            y -= x;
            k += 1;
        }
        a[k] = x + y;
        for (i, val) in a[..=k].iter().rev().enumerate() {
            results[[results_idx, i]] = *val as f32;
        }

        results_idx = match results_idx.checked_sub(1) {
            Some(val) => val,
            None => break,
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partitions_n1() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[[1.]]);
        assert_eq!(arr, partitions(1));
    }

    #[test]
    fn partitions_n3() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[[3., 0., 0.], [2., 1., 0.], [1., 1., 1.]]);
        assert_eq!(arr, partitions(3));
    }

    #[test]
    fn partitions_n5() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[
            [5., 0., 0., 0., 0.],
            [3., 2., 0., 0., 0.],
            [4., 1., 0., 0., 0.],
            [2., 2., 1., 0., 0.],
            [3., 1., 1., 0., 0.],
            [2., 1., 1., 1., 0.],
            [1., 1., 1., 1., 1.],
        ]);
        assert_eq!(arr, partitions(5));
    }

    #[test]
    fn partitions_n6() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[
            [6., 0., 0., 0., 0., 0.],
            [3., 3., 0., 0., 0., 0.],
            [4., 2., 0., 0., 0., 0.],
            [2., 2., 2., 0., 0., 0.],
            [5., 1., 0., 0., 0., 0.],
            [3., 2., 1., 0., 0., 0.],
            [4., 1., 1., 0., 0., 0.],
            [2., 2., 1., 1., 0., 0.],
            [3., 1., 1., 1., 0., 0.],
            [2., 1., 1., 1., 1., 0.],
            [1., 1., 1., 1., 1., 1.],
        ]);
        assert_eq!(arr, partitions(6));
    }

    #[test]
    fn partition_sums() {
        for (i, val) in N_PARTITIONS.iter().enumerate() {
            let parts = partitions(i);
            assert_eq!(*val, parts.shape()[0]);
            let rowsums = parts.sum_axis(Axis(1));
            for j in rowsums.iter() {
                assert_eq!(i as f32, *j);
            }
        }
    }

    #[test]
    fn most_likely_gt() {
        let partitions = arr2(&[[3., 0., 0.], [2., 1., 0.], [1., 1., 1.]]);
        let n_reads = 30;
        let cn = 3;
        let counts = arr1(&[20., 10., 0.]);

        assert_eq!(
            1,
            most_likely_partition_idx(&partitions, n_reads, cn, &counts).unwrap()
        );
    }

    #[test]
    fn multiple_most_likely_gts() {
        let partitions = arr2(&[[3., 0., 0.], [1., 1., 1.]]);
        let n_reads = 30;
        let cn = 3;
        let counts = arr1(&[20., 10., 0.]);

        let res = most_likely_partition_idx(&partitions, n_reads, cn, &counts);
        assert!(res.is_err());
    }

    #[test]
    fn valid_partitions() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[[3., 0., 0.], [2., 1., 0.], [1., 1., 1.]]);
        let counts: Array<f32, Dim<[usize; 1]>> = arr1(&[10., 8., 2.]);
        let threshold = 3.;
        let valid_partitions = find_valid_partition_idxs(&arr, &counts, threshold);

        assert_eq!(valid_partitions, vec![0, 1])
    }
}
