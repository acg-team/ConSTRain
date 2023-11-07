use std::collections::HashMap;

#[derive(Debug)]
pub struct TandemRepeat {
    pub reference_info: RepeatReferenceInfo,
    pub copy_number: u64,
    pub is_genotyped: bool,
    pub allele_lengths: Option<HashMap<i64, usize>>,
}

impl TandemRepeat {    
    pub fn set_genotyped(&mut self) {
        self.is_genotyped = true
    }
    pub fn set_allele_lengths(&mut self, new_allele_lengths: Option<HashMap<i64, usize>>) {
        self.allele_lengths = new_allele_lengths
    }
    pub fn set_cn(&mut self, new_cn: u64) {
        self.copy_number = new_cn
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