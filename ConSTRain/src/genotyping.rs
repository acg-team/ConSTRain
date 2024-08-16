//! # Estimating genotypes from allele length distributions
use anyhow::{bail, Context, Result};
use log::info;
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
    partitions_map: Arc<PartitionMap>,
) -> Result<()> {
    if tr_region.copy_number == 0 {
        bail!("cannot estimate genotype for locus with copy number 0");
    }
    let Some(n_mapped_reads) = tr_region.get_n_mapped_reads() else {
        bail!("cannot estimate genotype for locus with no reads mapped");
    };
    if n_mapped_reads < min_reads_per_allele * tr_region.copy_number {
        bail!("not enough reads were mapped to locus");
    }

    let (allele_lengths, mut obs_distr) = tr_region.allele_freqs_as_ndarrays(Some("freq"));

    let mut carry: Option<f32> = None;
    match obs_distr.len().cmp(&tr_region.copy_number) {
        Ordering::Less => obs_distr = utils::zero_pad_if_shorter(obs_distr, tr_region.copy_number),
        Ordering::Equal => (),
        Ordering::Greater => {
            carry = Some(obs_distr[tr_region.copy_number]); // keep track of this to check for edge cases later
            obs_distr = obs_distr.slice_move(s![..tr_region.copy_number])
        }
    };

    let possible_gts = partitions_map
        .get(&tr_region.copy_number)
        .with_context(|| {
            format!(
                "locus has copy number {}, for which no partitions exist",
                tr_region.copy_number
            )
        })?;

    let argmin = most_likely_gt_idx(
        &possible_gts,
        n_mapped_reads,
        tr_region.copy_number,
        &obs_distr,
    )?;
    let most_likely_gt = possible_gts.slice(s![argmin, ..]);

    let mut genotype: Vec<(i64, f32)> = allele_lengths
        .iter()
        .zip(most_likely_gt.iter())
        .filter_map(|(x, y)| if *y > 0. { Some((*x, *y)) } else { None })
        .collect();

    if let Some(val) = carry {
        if genotype.len() == tr_region.copy_number && obs_distr[tr_region.copy_number - 1] == val {
            bail!("multiple allele distributions are equally likely");
        }
    }

    genotype.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    tr_region.genotype = Some(genotype);

    Ok(())
}

/// Return the indexes of the most likely genotype to underly the observed allele distribution
fn most_likely_gt_idx(
    possible_gts: &Array<f32, Dim<[usize; 2]>>,
    n_mapped_reads: usize,
    copy_number: usize,
    obs_distr: &Array<f32, Dim<[usize; 1]>>,
) -> Result<usize> {
    // we assume each allele contributes the same number of reads to the observed
    // allele distribution. Thus, we find the expected number of reads per allele
    // by dividing the total number of mapped reads by the copy number
    let e_reads_per_allele = n_mapped_reads as f32 / copy_number as f32;
    let mut errors = possible_gts.mapv(|a| a * e_reads_per_allele);

    // Manhattan distance between observed and expected allele distributions
    errors = errors - obs_distr.slice(s![..copy_number]);
    errors.mapv_inplace(|x| x.abs());
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
            bail!("no most likely allele distribution was found. This should not happen")
        }
        Ordering::Equal => Ok(argmin[0]),
        Ordering::Greater => bail!("multiple allele distributions are equally likely"),
    }
}

pub type PartitionMap = HashMap<usize, Array<f32, Dim<[usize; 2]>>>;

pub fn make_partitions_map(copy_numbers: &[usize]) -> PartitionMap {
    info!("Generating possbile genotypes for copy number(s) {copy_numbers:?}");
    let mut map: PartitionMap = HashMap::new();
    for cn in copy_numbers {
        map.insert(*cn, partitions(*cn));
    }

    map
}

/// Generate all integer partitions of for n
///
/// Background: ([Wikipedia](https://en.wikipedia.org/wiki/Integer_partition), [Mathworld](https://mathworld.wolfram.com/Partition.html)).
/// All partitions are padded with zeroes to be of size 'n', and values are represented by f32s (because the partitions will be used to calculate errors later.)
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
        // Add the partitions in reversed order since we want
        // descending partitions (just because I like them better)
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
    use crate::repeat::RepeatReferenceInfo;

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
            most_likely_gt_idx(&partitions, n_reads, cn, &counts).unwrap()
        );
    }

    #[test]
    fn multiple_most_likely_gts() {
        let partitions = arr2(&[[2., 0.], [1., 1.]]);
        let cn = 2;
        let counts = arr1(&[3., 1.]);
        let n_reads = 4;
        
        assert!(most_likely_gt_idx(&partitions, n_reads, cn, &counts).is_err());
    }

    #[test]
    fn estimate_gt_pass1() {
        let ref_info = RepeatReferenceInfo {
            seqname: String::from("chr1"),
            start: 10,
            end: 31,
            period: 3,
            unit: String::from("CAG"),
        };
        let allele_lengths = HashMap::from([(10, 12.), (12, 10.), (13, 9.), (14, 3.)]);
        let mut tr = TandemRepeat {
            reference_info: ref_info,
            copy_number: 2,
            allele_lengths: Some(allele_lengths),
            genotype: None,
            skip: false,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));
        let target_gt: Vec<(i64, f32)> = vec![(10, 1.), (12, 1.)];

        assert!(estimate_genotype(&mut tr, 0, Arc::clone(&pmap)).is_ok());
        assert_eq!(tr.genotype, Some(target_gt));
    }

    #[test]
    fn estimate_gt_pass2() {
        let ref_info = RepeatReferenceInfo {
            seqname: String::from("chr1"),
            start: 10,
            end: 31,
            period: 3,
            unit: String::from("CAG"),
        };
        let allele_lengths = HashMap::from([(10, 20.), (12, 4.), (13, 4.), (14, 3.)]);
        let mut tr = TandemRepeat {
            reference_info: ref_info,
            copy_number: 2,
            allele_lengths: Some(allele_lengths),
            genotype: None,
            skip: false,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));
        let target_gt: Vec<(i64, f32)> = vec![(10, 2.)];

        assert!(estimate_genotype(&mut tr, 0, Arc::clone(&pmap)).is_ok());
        assert_eq!(tr.genotype, Some(target_gt));
    }

    #[test]
    fn estimate_gt_skip() {
        let ref_info = RepeatReferenceInfo {
            seqname: String::from("chr1"),
            start: 10,
            end: 31,
            period: 3,
            unit: String::from("CAG"),
        };
        let allele_lengths = HashMap::from([(10, 12.), (12, 10.), (13, 10.), (14, 3.)]);
        let mut tr = TandemRepeat {
            reference_info: ref_info,
            copy_number: 2,
            allele_lengths: Some(allele_lengths),
            genotype: None,
            skip: false,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));

        assert!(estimate_genotype(&mut tr, 0, Arc::clone(&pmap)).is_err());
        assert_eq!(tr.genotype, None);
    }
}
