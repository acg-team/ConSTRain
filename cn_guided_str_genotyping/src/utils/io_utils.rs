use std::{fs::File, io::BufReader};

use csv::ReaderBuilder;
use serde_json::Value;

use crate::repeat::{TandemRepeat, RepeatReferenceInfo};

pub fn trs_from_bed(
    bed_path    : &String, 
    ploidy_path : &String,
    tr_regions  : &mut Vec<TandemRepeat>
) {
    let ploidy = load_ploidy(ploidy_path);
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .flexible(true)
        .from_path(bed_path).unwrap();    

    for result in bed_reader.deserialize() {
        let record: RepeatReferenceInfo = result.unwrap();
        let mut tr_region = TandemRepeat{
            reference_info  : record, 
            copy_number     : 0, 
            is_genotyped    : false, 
            allele_lengths  : None
        };

        let cn = ploidy[tr_region.reference_info.seqname.as_str()].as_u64().unwrap();
        tr_region.set_cn(cn);

        tr_regions.push(tr_region);
    }
}

fn load_ploidy(path: &String) -> Value {
    let f = File::open(path).unwrap();
    let reader = BufReader::new(f);
    let v: Value = serde_json::from_reader(reader).unwrap();

    v
}