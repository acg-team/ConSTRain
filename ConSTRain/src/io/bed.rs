use std::collections::{HashMap, HashSet};

use anyhow::{bail, Context, Result};
use csv::ReaderBuilder;
use log::{error, info, warn};

use crate::{
    io::json::ploidy_from_json,
    repeat::{RepeatReferenceInfo, TandemRepeat},
    utils::CopyNumberVariant,
};

/// Read tandem repeats from `repeats` (a bed file) and set their baseline copy number based
/// on the genome architecture specified in `ploidy` (a json). Optionally, parse copy number
/// variants from `cnvs` (a bed file) and update the copy number of tandem repeats located
/// in those regions.
pub fn parse_tandem_repeats(
    repeats: &str,
    ploidy: &str,
    cnvs: Option<&str>,
) -> Result<(Vec<TandemRepeat>, Vec<usize>)> {
    let mut observed_copy_numbers: HashSet<usize> = HashSet::new();
    let mut tr_regions: Vec<TandemRepeat> = Vec::new();

    read_trs(repeats, ploidy, &mut tr_regions, &mut observed_copy_numbers)?;
    info!("Read {} TR regions from {}", tr_regions.len(), repeats);

    if let Some(cnv_file) = cnvs {
        let cnv_map = read_cnvs(cnv_file, &mut observed_copy_numbers)?;
        info!("Read CNVs from {}", cnv_file);
        for tr_region in &mut tr_regions {
            if let Some(cnv_vec) = cnv_map.get(&tr_region.reference_info.seqname) {
                if let Err(e) = tr_region.set_cn_from_cnvs(cnv_vec) {
                    error!(
                        "Error setting copy number from CNV for locus {}: {e:?}",
                        tr_region.reference_info.get_fetch_definition_s()
                    );
                    tr_region.set_cn(0); // set CN to 0. Or should we continue with baseline CN here?
                    continue;
                };
            };
        }
        info!("Updated TR copy numbers using CNVs");
    }

    let mut observed_copy_numbers: Vec<usize> = observed_copy_numbers
        .iter()
        .filter_map(|x| if 0 < *x && *x <= 20 { Some(*x) } else { None }) // ConSTRain supports copy numbers from 1 to 20 for now
        .collect();
    observed_copy_numbers.sort();

    Ok((tr_regions, observed_copy_numbers))
}

/// Read tandem repeat regions specified in the bed file at `bed_path` into `tr_buffer`.
/// Copy numbers of tandem repeats are set based on the values specified in `ploidy_path`, which should be a
/// json file mapping contig identifiers to copy number values. All copy number values that are observed while
/// reading TR regions are added to `cn_buffer`, which is used later to generate partitions for only the
/// relevant (i.e., observed) copy numbers.
fn read_trs(
    bed_path: &str,
    ploidy_path: &str,
    tr_buffer: &mut Vec<TandemRepeat>,
    observed_cn_buffer: &mut HashSet<usize>,
) -> Result<()> {
    let ploidy = ploidy_from_json(ploidy_path)?;
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(bed_path)
        .with_context(|| format!("Could not read bed file {bed_path}"))?;

    for result in bed_reader.deserialize() {
        let ref_info: RepeatReferenceInfo =
            result.with_context(|| format!("Failed to deserialize bed record in {bed_path}"))?;

        // Get base copy number of contig from ploidy json.
        // If the contig does not exist in the json, default to CN 0 (which means TRs on this contig won't be genotyped)
        let cn = ploidy[&ref_info.seqname].as_u64().unwrap_or(0) as usize;
        if cn == 0 {
            warn!(
                "Contig for TR '{}' was not found in the ploidy json file. Setting CN to 0",
                ref_info.get_fetch_definition_s()
            );
        }

        let tr_region = TandemRepeat {
            reference_info: ref_info,
            copy_number: cn,
            allele_lengths: None,
            genotype: None,
        };
        observed_cn_buffer.insert(cn);
        tr_buffer.push(tr_region);
    }

    Ok(())
}

/// Read copy number variants specified in the bed file at `cnv_path` into a HashMap where keys are contig
/// names and values are vectors of CNVs. Ensure that CNVs are coordinate sorted.
/// All copy number values that are observed while reading TR regions are added to `cn_buffer`,
/// which is used later to generate partitions for only the relevant (i.e., observed) copy numbers.
fn read_cnvs(
    cnv_path: &str,
    observed_cn_buffer: &mut HashSet<usize>,
) -> Result<HashMap<String, Vec<CopyNumberVariant>>> {
    let mut bed_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(cnv_path)
        .with_context(|| format!("Could not read bed file {cnv_path}"))?;

    let mut cnv_map: HashMap<String, Vec<CopyNumberVariant>> = HashMap::new();
    for result in bed_reader.deserialize() {
        let cnv: CopyNumberVariant =
            result.with_context(|| format!("Failed to deserialize bed record in {cnv_path}"))?;
        observed_cn_buffer.insert(cnv.cn);
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

    Ok(cnv_map)
}
