/// Struct to represent Copy Number Variants
#[derive(Debug, serde::Deserialize)]
pub struct CopyNumberVariant {
    pub seqname: String,
    // CNV struct follows 0-based half-open coordinate system: [start, end)
    pub start: i64,
    pub end: i64,
    pub cn: usize,
}
