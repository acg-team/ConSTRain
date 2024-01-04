use crate::repeat::TandemRepeat;
use crate::utils::N_PARTITIONS;

use ndarray::prelude::*;
use std::{collections::HashMap, sync::Arc};

pub fn estimate_genotype(
    tr_region: &mut TandemRepeat,
    min_reads_per_allele: usize,
    partitions_map: Arc<HashMap<usize, Array<f32, Dim<[usize; 2]>>>>,
) {
    let n_mapped_reads = match tr_region.get_n_mapped_reads() {
        Some(n) => n,
        None => return, // No reads were mapped to this locus, refuse to estimate
    };
    if n_mapped_reads < min_reads_per_allele * tr_region.copy_number {
        // Not enough reads were mapped to this locus, refuse to estimate
        return;
    }

    let (allele_lengths, mut counts) = tr_region.allele_counts_as_ndarrays();
    // might need this error_constant at some point if trying to infer the CN from read distribution
    // let mut error_constant = 0.;
    if counts.len() < tr_region.copy_number {
        counts = zero_pad_if_shorter(
            counts,
            tr_region.copy_number,
        );
    } else if counts.len() > tr_region.copy_number {
        // error_constant += counts.slice(s![tr_region.copy_number..]).sum();
        counts = counts.slice_move(s![..tr_region.copy_number]);
    }
    

    let partitions = match partitions_map.get(&tr_region.copy_number) {
        Some(comp) => comp,
        None => return,
    };

    let threshold_val = n_mapped_reads as f32 / tr_region.copy_number as f32 * 0.5;
    let valid_partition_idxs = find_valid_partition_idxs(&partitions, &counts, threshold_val);

    if valid_partition_idxs.len() == 0 {
        // Not a single allele length observed with frequency over threshold, refuse to estimate
        return;        
    } else if valid_partition_idxs.len() == 1 {
        // Only one allele length observed with frequency over threshold, we can exit early
        let genotype: Vec<(i64, f32)> = vec![(allele_lengths[0], tr_region.copy_number as f32)];
        tr_region.genotype = Some(genotype);
        return; 
    }

    let valid_partitions = partitions.select(Axis(0), &valid_partition_idxs);

    let argmin = most_likely_allele_distribution(
        &valid_partitions,
        n_mapped_reads,
        tr_region.copy_number,
        &counts,
    );
    let allele_distribution = valid_partitions.slice(s![argmin, ..]);

    let genotype: Vec<(i64, f32)> = allele_lengths
        .iter()
        .zip(allele_distribution.iter())
        .filter_map(|(x, y)| if *y > 0. { Some((*x, *y)) } else { None })
        .collect();
    tr_region.genotype = Some(genotype);
}

fn zero_pad_if_shorter(
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

fn find_valid_partition_idxs(
    partitions: &Array<f32, Dim<[usize; 2]>>,
    counts: &Array<f32, Dim<[usize; 1]>>,
    threshold_value: f32,
) -> Vec<usize> {
    let mut threshold_idx = 0;
    for count in counts.iter() {
        if *count > threshold_value {
            threshold_idx += 1;
        } else {
            break;
        }
    }

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

fn most_likely_allele_distribution(
    partitions: &Array<f32, Dim<[usize; 2]>>,
    n_mapped_reads: usize,
    copy_number: usize,
    counts: &Array<f32, Dim<[usize; 1]>>,
) -> usize {
    let mut errors = partitions.mapv(|a| a * (n_mapped_reads / copy_number) as f32);
    errors = errors - counts.slice(s![..copy_number]);
    errors.mapv_inplace(|x| x.powi(2));

    let error_sums = errors.sum_axis(Axis(1));
    let argmin = error_sums
        .iter()
        .enumerate()
        .min_by(|x, y| x.1.partial_cmp(y.1).unwrap())
        .unwrap();

    argmin.0
}

pub fn partitions(n: usize) -> Array<f32, Dim<[usize; 2]>> {
    // https://jeromekelleher.net/category/combinatorics
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
            k += 1
        }
        a[k] = x + y;
        for (i, val) in a[..k + 1].iter().rev().enumerate() {
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
    fn zero_padding() {
        let left: Array<f32, Dim<[usize; 1]>> = arr1(&[13., 12., 11.]);
        let right: Array<f32, Dim<[usize; 1]>> = arr1(&[13., 12., 11.]);
        assert_eq!(left, zero_pad_if_shorter(right, 3));

        let left: Array<f32, Dim<[usize; 1]>> = arr1(&[13., 12., 11., 0., 0.]);
        let right: Array<f32, Dim<[usize; 1]>> = arr1(&[13., 12., 11.]);
        assert_eq!(left, zero_pad_if_shorter(right, 5));
    }

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
    fn partitions_test() {
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
        let partitions: Array<f32, Dim<[usize; 2]>> =
            arr2(&[[3., 0., 0.], [2., 1., 0.], [1., 1., 1.]]);
        let n_reads = 30;
        let cn = 3;
        let counts: Array<f32, Dim<[usize; 1]>> = arr1(&[20., 10., 0.]);

        assert_eq!(
            1,
            most_likely_allele_distribution(&partitions, n_reads, cn, &counts)
        );
    }
}
