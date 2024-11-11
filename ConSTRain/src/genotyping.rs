//! # Estimating genotypes from allele length distributions
use anyhow::{bail, Context, Result};
use log::info;
use ndarray::{prelude::*, Slice};
use std::{cmp::Ordering, collections::HashMap, sync::Arc};

use crate::{
    repeat::TandemRepeat,
    utils::{self, VcfFilter, N_PARTITIONS},
};

/// Estimate the most likely underlying genotype that produced the
/// observed allele distribution stored on the `tr_region` struct.
pub fn estimate_genotype(
    tr_region: &mut TandemRepeat,
    min_norm_depth: f32,
    max_norm_depth: Option<f32>,
    partitions_map: Arc<PartitionMap>,
) -> Result<()> {
    if tr_region.copy_number == 0 {
        tr_region.filter = VcfFilter::CnZero;
        bail!("cannot estimate genotype for locus with copy number 0");
    }
    let Some(norm_depth) = tr_region.get_norm_depth() else {
        tr_region.filter = VcfFilter::DpZero;
        bail!("no reads were mapped to locus");
    };
    if norm_depth < min_norm_depth {
        tr_region.filter = VcfFilter::DpOor;
        bail!("locus normalised depth is too low");
    }
    if let Some(max_norm_depth) = max_norm_depth {
        if norm_depth > max_norm_depth {
            tr_region.filter = VcfFilter::DpOor;
            bail!("locus normalised depth is too high");
        }
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

    let Some(possible_gts) = partitions_map.get(&tr_region.copy_number) else {
        tr_region.filter = VcfFilter::CnOor;
        bail!(
            "locus has copy number {}, for which no partitions were generated",
            tr_region.copy_number
        )
    };

    let min = match most_likely_gt_idx(
        possible_gts,
        tr_region
            .get_n_mapped_reads()
            .context("no reads were mapped to locus (should NOT be possible here)")?,
        tr_region.copy_number,
        &obs_distr,
    ) {
        Ok(min) => min,
        Err(e) => {
            tr_region.filter = VcfFilter::AmbGt;
            return Err(e);
        }
    };

    let most_likely_gt = possible_gts.slice(s![min, ..]);
    // If there are any 'plateaus' in observed distribution, check that they all have the same allele frequency
    if let Some(plateaus) = find_plateaus(&obs_distr) {
        for p in plateaus {
            let subarr = most_likely_gt.slice(s![p]);
            if !subarr.slice(s![1..]).iter().all(|v| *v == subarr[0]) {
                tr_region.filter = VcfFilter::AmbGt;
                bail!("multiple allele distributions are equally likely");
            }
        }
    }
    let mut genotype: Vec<(i64, f32)> = allele_lengths
        .iter()
        .zip(most_likely_gt.iter())
        .filter_map(|(x, y)| if *y > 0. { Some((*x, *y)) } else { None })
        .collect();

    if let Some(val) = carry {
        if genotype.len() == tr_region.copy_number && obs_distr[tr_region.copy_number - 1] == val {
            tr_region.filter = VcfFilter::AmbGt;
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

    let min: Vec<usize> = error_sums
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

    match min.len().cmp(&1) {
        Ordering::Less => {
            bail!("no most likely allele distribution was found. This should not happen")
        }
        Ordering::Equal => Ok(min[0]),
        Ordering::Greater => bail!("multiple allele distributions are equally likely"),
    }
}

/// Check if there are allele lengths in the observed allele length distribution
/// that are represented by exactly the same number of reads. We call these 'plateaus'.
/// Under ConSTRain's assumptions, all alleles in a plateau should have the
/// same allele frequency in the final genotype.
fn find_plateaus(arr: &Array<f32, Dim<[usize; 1]>>) -> Option<Vec<Slice>> {
    let mut plateaus: Vec<Slice> = Vec::new();
    let mut current: Option<(usize, usize)> = None;

    for (i, val) in arr.slice(s![1..]).iter().enumerate() {
        if arr[i] == *val {
            if let Some(ref mut p) = current {
                // extend the current plateau
                p.1 += 1;
            } else {
                current = Some((i, i + 2));
            }
        } else {
            if let Some(p) = current {
                // reached end of plateau, push and reset
                plateaus.push(Slice::new(p.0 as isize, Some(p.1 as isize), 1));
                current = None;
            }
        }
    }

    if let Some(p) = current {
        // if there was still a `current` at this point, we can slice to the end
        plateaus.push(Slice::new(p.0 as isize, None, 1));
    }

    if plateaus.len() > 0 {
        return Some(plateaus);
    }

    None
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
    fn plateaus() {
        let arr = Array::from_vec(vec![48., 32., 32., 8., 4., 4., 4., 2., 2.]);
        let expect = Some(vec![
            Slice {
                start: 1,
                end: Some(3),
                step: 1,
            },
            Slice {
                start: 4,
                end: Some(7),
                step: 1,
            },
            Slice {
                start: 7,
                end: None,
                step: 1,
            },
        ]);

        assert_eq!(find_plateaus(&arr), expect);
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
    fn multiple_most_likely_gts2() {
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
            filter: utils::VcfFilter::Pass,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));

        assert!(estimate_genotype(&mut tr, 0., None, Arc::clone(&pmap)).is_err());
        assert_eq!(tr.genotype, None);
    }

    #[test]
    fn multiple_most_likely_gts3() {
        // {12: 4, 13: 32, 14: 4}
        // chr1_155186468
        let ref_info = RepeatReferenceInfo {
            seqname: String::from("chr1"),
            start: 155186468,
            end: 155186492,
            period: 2,
            unit: String::from("GT"),
        };

        let allele_lengths = HashMap::from([(12, 4.), (13, 32.), (14, 4.)]);
        let mut tr = TandemRepeat {
            reference_info: ref_info,
            copy_number: 4,
            allele_lengths: Some(allele_lengths),
            genotype: None,
            filter: utils::VcfFilter::Pass,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));

        assert!(estimate_genotype(&mut tr, 0., None, Arc::clone(&pmap)).is_err());
        assert_eq!(tr.genotype, None);
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
            filter: utils::VcfFilter::Pass,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));
        let target_gt: Vec<(i64, f32)> = vec![(10, 1.), (12, 1.)];

        assert!(estimate_genotype(&mut tr, 0., None, Arc::clone(&pmap)).is_ok());
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
            filter: utils::VcfFilter::Pass,
        };
        let pmap = Arc::new(make_partitions_map(&vec![tr.copy_number]));
        let target_gt: Vec<(i64, f32)> = vec![(10, 2.)];

        assert!(estimate_genotype(&mut tr, 0., None, Arc::clone(&pmap)).is_ok());
        assert_eq!(tr.genotype, Some(target_gt));
    }
}
