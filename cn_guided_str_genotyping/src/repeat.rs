use ndarray::prelude::*;
use std::collections::HashMap;

#[derive(Debug)]
pub struct TandemRepeat {
    pub reference_info: RepeatReferenceInfo,
    pub copy_number: u64,
    pub allele_lengths: Option<HashMap<i64, f32>>,
    pub n_mapped_reads: u64,
    pub genotype: Option<Vec<(i64, f32)>>,
}

impl TandemRepeat {
    pub fn has_coverage(&self) -> bool {
        self.allele_lengths.is_some()
    }
    pub fn set_allele_lengths(&mut self, new_allele_lengths: Option<HashMap<i64, f32>>) {
        self.allele_lengths = new_allele_lengths
    }
    pub fn set_cn(&mut self, new_cn: u64) {
        self.copy_number = new_cn
    }
    pub fn allele_lengths_as_ndarrays(&self) -> (Array::<i64, Dim<[usize; 1]>>, Array::<f32, Dim<[usize; 1]>>) {
        let mut count_vec: Vec<(&i64, &f32)> =
            self.allele_lengths.as_ref().unwrap().iter().collect();
        count_vec.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());

        let allele_lengths = Array::<i64, _>::from_vec(count_vec.iter().map(|a| *a.0).collect());
        let counts = Array::<f32, _>::from_vec(count_vec.iter().map(|a| *a.1).collect());

        return(allele_lengths, counts)
    }
}

#[derive(Debug, serde::Deserialize)]
pub struct RepeatReferenceInfo {
    pub seqname: String,
    // Tandem repeat struct follows 0-based half-open coordinate system: [start, end)
    pub start: i64,
    pub end: i64,
    pub period: i64,
    pub unit: String,
}

impl RepeatReferenceInfo {
    pub fn get_fetch_definition(&self) -> (&str, i64, i64) {
        (self.seqname.as_str(), self.start, self.end)
    }
}
