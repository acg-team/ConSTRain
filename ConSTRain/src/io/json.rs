use std::{fs::File, io::BufReader};

use anyhow::{Context, Result};
use log::debug;
use serde_json::Value;

/// Read contig-level baseline copy number values from a json file at `path`.
/// The json should contain contig names as keys and integer copy numbers as values, e.g.:
/// `
/// {
///     "chr1": 2,
///     ... other chromosomes ...
///     "chrY": 0
/// }
/// `
pub fn ploidy_from_json(path: &str) -> Result<Value> {
    let file = File::open(path).with_context(|| format!("Could not read json {path}"))?;
    let reader = BufReader::new(file);
    let result = serde_json::from_reader(reader)
        .with_context(|| format!("Could not deserialize json {path}"))?;
    Ok(result)
}

pub fn get_cn(ploidy: &Value, contig: &str) -> Option<usize> {
    if let Some(val) = ploidy[contig].as_u64() {
        Some(val as usize)
    } else {
        debug!("Contig for '{contig}' was not found in the ploidy json file.",);
        None
    }
}
