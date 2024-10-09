//! # CIGAR Utils
//!
//! Functions to check which category CIGAR operations fall into,
//! which informs how they should be interpreted during genotyping.
use rust_htslib::bam::record::Cigar;

/// Check if the provided cigar operation `cigar` advances the
/// position in the reference sequence.
pub fn consumes_ref(cigar: &Cigar) -> bool {
    matches!(
        cigar,
        Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_)
    )
}

/// Check if the provided cigar operation `cigar` advances the
/// position in the query sequence.
pub fn consumes_query(cigar: &Cigar) -> bool {
    matches!(
        cigar,
        Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_)
    )
}

/// Used when in a TR region to see if the provided cigar operation
/// `cigar` advances the TR length.
/// (same as [`consumes_query`] but excluding the [`Cigar::SoftClip`] variant)
pub fn advances_tr_len(cigar: &Cigar) -> bool {
    matches!(
        cigar,
        Cigar::Match(_) | Cigar::Ins(_) | Cigar::Equal(_) | Cigar::Diff(_)
    )
}
