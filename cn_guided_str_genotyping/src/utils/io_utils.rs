use csv::ReaderBuilder;

use crate::utils::ReferenceRepeat;

pub fn trs_from_bed(
    bed_path    : &String, 
    tr_regions : &mut Vec<ReferenceRepeat>
) {
    let mut reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .flexible(true)
        .from_path(bed_path).unwrap();

    for result in reader.deserialize() {
        let mut record: ReferenceRepeat = result.unwrap();
        // record.start -= 1;
        // record.end -= 1;
        tr_regions.push(record);
    }
}