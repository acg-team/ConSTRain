use anyhow::Result;
use log::debug;
use serde_json::Value;

use crate::io::json::read_karyotype;

pub struct Karyotype {
    ploidies: Value,
}

impl Karyotype {
    pub fn from_json(path: &str) -> Result<Karyotype> {
        Ok(Karyotype {
            ploidies: read_karyotype(path)?,
        })
    }
    pub fn get_ploidy(&self, contig: &str) -> Option<usize> {
        if let Some(val) = self.ploidies[contig].as_u64() {
            Some(val as usize)
        } else {
            debug!("Contig '{contig}' was not present in the karyotype",);
            None
        }
    }
}
