use crate::repeat::TandemRepeat;

use ndarray::prelude::*;

pub fn estimate_genotype(tr_region: &mut TandemRepeat, min_reads_per_allele: u64) {
    if !tr_region.has_coverage() 
        || tr_region.n_mapped_reads <= min_reads_per_allele * tr_region.copy_number {
        return;
    }

    let (allele_lengths, mut counts) = tr_region.allele_lengths_as_ndarrays();
    counts = zero_pad_if_shorter(counts, tr_region.copy_number as usize);

    let compositions = descending_weak_compositions(tr_region.copy_number as usize);    
    let mut errors = compositions.mapv(|a| a * (tr_region.n_mapped_reads / tr_region.copy_number) as f32);
    errors = errors - counts.slice(s![..tr_region.copy_number as usize]);
    errors.mapv_inplace(|x| x.powi(2));    

    let error_sums = errors.sum_axis(Axis(1));
    let min = error_sums.iter().enumerate().min_by(|x, y| x.1.partial_cmp(y.1).unwrap()).unwrap();
    let allele_distribution = compositions.slice(s![min.0,..]);

    let genotype: Vec<(i64, f32)> = allele_lengths
                                .iter()
                                .zip(allele_distribution.iter())
                                .filter_map(|x| if *x.1 > 0. { Some((x.0.to_owned(), x.1.to_owned() )) } else { None })
                                .collect();
    tr_region.genotype = Some(genotype);
}

pub fn zero_pad_if_shorter(a: Array::<f32, Dim<[usize; 1]>>, min_len: usize) -> Array::<f32, Dim<[usize; 1]>> {
    if a.len() >= min_len {
        return a;
    }
    let mut padded_a = Array::<f32, _>::zeros(min_len);
    padded_a.slice_mut(s![..a.len()]).assign(&a);
    padded_a
}

pub fn descending_weak_compositions(n: usize) -> Array<f32, Dim<[usize; 2]>> {
    let mut results: Vec<f32> = Vec::new();
    let mut composition = Array::<usize, _>::zeros(n);
    composition[0] = n;

    results.extend_from_slice(&composition.mapv(|x| x as f32).to_vec());

    let mut idx = 0;
    let mut n_results = 1;
    loop {
        if composition[idx] > 1 {
            composition[idx] -= 1;
            let sum_right = composition.slice(s![idx + 1 ..]).sum() + 1; // include the count that we just decremented from composition[idx]
            if sum_right >= 2 {
                let current_val = composition[idx];
                let divmod = (sum_right / current_val, sum_right % current_val);
                
                // handle quotient (could be optimized more)
                composition.slice_mut(s![idx + 1 ..]).fill(0);
                composition.slice_mut(s![idx + 1 .. idx + 1 + divmod.0]).fill(current_val);
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