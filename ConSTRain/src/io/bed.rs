use std::collections::{HashMap, HashSet};

use anyhow::{bail, Context, Result};
use csv::ReaderBuilder;
use log::info;

use crate::{
    cnv::CopyNumberVariant,
    karyotype::Karyotype,
    repeat::{RepeatReferenceInfo, TandemRepeat},
};

/// Read tandem repeat regions specified in the bed file at `bed_path` into `tr_buffer`.
/// Copy numbers of tandem repeats are set based on the values specified in `ploidy_path`, which should be a
/// json file mapping contig identifiers to copy number values. All copy number values that are observed while
/// reading TR regions are added to `cn_buffer`, which is used later to generate partitions for only the
/// relevant (i.e., observed) copy numbers.
pub fn read_trs(
    bed_path: &str,
    karyotype: &Karyotype,
    tr_buffer: &mut Vec<TandemRepeat>,
    observed_cn_buffer: &mut HashSet<usize>,
) -> Result<()> {
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(bed_path)
        .with_context(|| format!("Could not read bed file {bed_path}"))?;

    for result in bed_reader.deserialize() {
        let ref_info: RepeatReferenceInfo =
            result.with_context(|| format!("Failed to deserialize bed record in {bed_path}"))?;

        let mut tr = TandemRepeat {
            reference_info: ref_info,
            copy_number: 0, // placeholder copy number value
            allele_lengths: None,
            genotype: None,
            skip: false,
        };

        tr.set_cn_from_karyotpe(karyotype);
        observed_cn_buffer.insert(tr.copy_number);
        tr_buffer.push(tr);
    }

    info!("Read {} TR regions from {bed_path}", tr_buffer.len());
    Ok(())
}

/// Read copy number variants specified in the bed file at `cnv_path` into a HashMap where keys are contig
/// names and values are vectors of CNVs. Ensure that CNVs are coordinate sorted.
/// All copy number values that are observed while reading TR regions are added to `cn_buffer`,
/// which is used later to generate partitions for only the relevant (i.e., observed) copy numbers.
pub fn read_cnvs(
    cnv_path: &str,
    observed_cn_buffer: &mut HashSet<usize>,
) -> Result<HashMap<String, Vec<CopyNumberVariant>>> {
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(cnv_path)
        .with_context(|| format!("Could not read bed file {cnv_path}"))?;

    let mut cnv_map: HashMap<String, Vec<CopyNumberVariant>> = HashMap::new();
    let mut n = 0;
    for result in bed_reader.deserialize() {
        let cnv: CopyNumberVariant =
            result.with_context(|| format!("Failed to deserialize bed record in {cnv_path}"))?;
        observed_cn_buffer.insert(cnv.cn);
        n += 1;
        if let Some(cnv_vec) = cnv_map.get_mut(&cnv.seqname) {
            let prev_cnv = &cnv_vec[cnv_vec.len() - 1];
            if prev_cnv.end > cnv.start && prev_cnv.seqname == cnv.seqname {
                bail!("CNV file {cnv_path} is not coordinate sorted");
            }
            cnv_vec.push(cnv);
        } else {
            cnv_map.insert(cnv.seqname.clone(), vec![cnv]);
        }
    }

    info!("Read {n} CNVS from {cnv_path}");

    Ok(cnv_map)
}
