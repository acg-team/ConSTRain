use crate::repeat::TandemRepeat;

use ndarray::prelude::*;
use std::{collections::HashMap, sync::Arc};

pub fn estimate_genotype(
    tr_region: &mut TandemRepeat,
    min_reads_per_allele: usize,
    compositions_map: Arc<HashMap<usize, Array<f32, Dim<[usize; 2]>>>>,
) {
    let n_mapped_reads = match tr_region.get_n_mapped_reads() {
        Some(n) => n,
        None => return,
    };
    if n_mapped_reads < min_reads_per_allele * tr_region.copy_number {
        return;
    }

    let (allele_lengths, mut counts) = tr_region.allele_counts_as_ndarrays();

    let threshold_count = n_mapped_reads as f32 / tr_region.copy_number as f32 * 0.5;

    // Determine idx of last position in counts that is higher than the threshold
    let mut threshold_idx = 0;
    for count in counts.iter() {
        if *count > threshold_count {
            threshold_idx += 1;
        } else {
            break;
        }
    }
    if threshold_idx == 0 {
        return;
    }
    counts = zero_pad_if_shorter(
        counts.slice_move(s![..threshold_idx]),
        tr_region.copy_number,
    );

    if allele_lengths.len() == 1 {
        // Only one allele length observed so we can exit early
        let genotype: Vec<(i64, f32)> = vec![(allele_lengths[0], tr_region.copy_number as f32)];
        tr_region.genotype = Some(genotype);
        return;
    }

    let compositions = match compositions_map.get(&tr_region.copy_number) {
        Some(comp) => comp,
        None => return,
    };

    // TODO: determine valid rows here using threshold_idx, pass to most_likely_allele_distribution
    // let valid_rows: Vec<usize> = compositions
    //     .slice(s![.., ..threshold_idx])
    //     .sum_axis(Axis(1))
    //     .iter()
    //     .enumerate()
    //     .filter_map(
    //         |(idx, val)| {
    //             if *val as usize == tr_region.copy_number {
    //                 Some(idx)
    //             } else {
    //                 None
    //             }
    //         },
    //     )
    //     .collect();
    // let valid_compositions = compositions.select(Axis(0), &valid_rows);

    let argmin = most_likely_allele_distribution(
        compositions,
        n_mapped_reads,
        tr_region.copy_number,
        &counts,
    );
    let allele_distribution = compositions.slice(s![argmin, ..]);

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

fn most_likely_allele_distribution(
    compositions: &Array<f32, Dim<[usize; 2]>>,
    n_mapped_reads: usize,
    copy_number: usize,
    counts: &Array<f32, Dim<[usize; 1]>>,
) -> usize {
    let mut errors = compositions.mapv(|a| a * (n_mapped_reads / copy_number) as f32);
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

pub fn descending_weak_compositions(n: usize) -> Array<f32, Dim<[usize; 2]>> {
    if n == 0 {
        return arr2(&[[]]);
    }
    let mut results: Vec<f32> = Vec::new();
    let mut composition = Array::<usize, _>::zeros(n);
    composition[0] = n;

    results.extend_from_slice(&composition.mapv(|x| x as f32).to_vec());

    let mut idx = 0;
    let mut n_results = 1;
    loop {
        if composition[idx] > 1 {
            composition[idx] -= 1;
            let sum_right = composition.slice(s![idx + 1..]).sum() + 1; // include the count that we just decremented from composition[idx]
            if sum_right >= 2 {
                let current_val = composition[idx];
                let divmod = (sum_right / current_val, sum_right % current_val);

                // handle quotient (could be optimized more)
                composition.slice_mut(s![idx + 1..]).fill(0);
                composition
                    .slice_mut(s![idx + 1..idx + 1 + divmod.0])
                    .fill(current_val);
                idx += divmod.0;

                // handle remainder (if there is one)
                if divmod.1 > 0 {
                    composition[idx + 1] = divmod.1;
                    idx += 1
                }
            } else {
                composition[idx + 1] += 1;
                idx += 1;
            }

            if composition[idx] <= composition[idx - 1] {
                // current position is lower than previous: we've
                // found a descending weak composition!
                results.extend_from_slice(&composition.mapv(|x| x as f32).to_vec());
                n_results += 1;
            }
        } else if idx > 0 {
            idx -= 1;
        } else {
            break;
        }
    }

    Array2::from_shape_vec((n_results, n), results).unwrap()
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
    fn compositions_n1() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[[1.]]);
        assert_eq!(arr, descending_weak_compositions(1));
    }

    #[test]
    fn compositions_n3() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[[3., 0., 0.], [2., 1., 0.], [1., 1., 1.]]);
        assert_eq!(arr, descending_weak_compositions(3));
    }

    #[test]
    fn compositions_n5() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[
            [5., 0., 0., 0., 0.],
            [4., 1., 0., 0., 0.],
            [3., 2., 0., 0., 0.],
            [3., 1., 1., 0., 0.],
            [2., 2., 1., 0., 0.],
            [2., 1., 1., 1., 0.],
            [1., 1., 1., 1., 1.],
        ]);
        assert_eq!(arr, descending_weak_compositions(5));
    }

    #[test]
    fn compositions_n6() {
        let arr: Array<f32, Dim<[usize; 2]>> = arr2(&[
            [6., 0., 0., 0., 0., 0.],
            [5., 1., 0., 0., 0., 0.],
            [4., 2., 0., 0., 0., 0.],
            [4., 1., 1., 0., 0., 0.],
            [3., 3., 0., 0., 0., 0.],
            [3., 2., 1., 0., 0., 0.],
            [3., 1., 1., 1., 0., 0.],
            [2., 2., 2., 0., 0., 0.],
            [2., 2., 1., 1., 0., 0.],
            [2., 1., 1., 1., 1., 0.],
            [1., 1., 1., 1., 1., 1.],
        ]);
        assert_eq!(arr, descending_weak_compositions(6));
    }

    #[test]
    fn most_likely_gt() {
        let compositions: Array<f32, Dim<[usize; 2]>> =
            arr2(&[[3., 0., 0.], [2., 1., 0.], [1., 1., 1.]]);
        let n_reads = 30;
        let cn = 3;
        let counts: Array<f32, Dim<[usize; 1]>> = arr1(&[20., 10., 0.]);

        assert_eq!(
            1,
            most_likely_allele_distribution(&compositions, n_reads, cn, &counts)
        );
    }
}
