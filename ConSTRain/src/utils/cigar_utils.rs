use rust_htslib::bam::record::Cigar;

pub fn cigar_consumes_ref(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
            true
        }
        _ => false,
    }
}

pub fn cigar_consumes_query(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
            true
        }
        _ => false,
    }
}

pub fn cigar_advances_tr_len(cigar: &Cigar) -> bool {
    match cigar {
        Cigar::Match(_)     |
        Cigar::Ins(_)       |
        // Cigar::SoftClip(_)  | // softclip consumes reference but shouldn't advance TR length (right?)
        Cigar::Equal(_)     |
        Cigar::Diff(_)      => true,
        _                   => false,
    }
}
