//! # Reimplementations of rust htslib private functions to use in ConSTRain
//!
//! Functions in this module are prefixed with 'rhtslib_'. They  are private functions in
//! rust_htslib and are copied from there to be used in this binary.
//! It would be nicer to use rust_htslib's structs (e.g., IndexedReader) for reading
//! alignment files, but those structs cause segfaults when reading CRAM files (not BAM, interestingly)
//! when they are dropped, even on trivial tests -> submit issue on GitHub?
use anyhow::{bail, Context, Result};
use rust_htslib::{
    bam::{self, Record},
    errors::Error as htsError,
    htslib::{self},
};
use std::{ffi, path::Path};

pub fn rhtslib_from_path<P: AsRef<Path>>(alignment_path: P) -> Result<*mut htslib::htsFile> {
    let alignment_path = rust_htslib::utils::path_as_bytes(alignment_path, true)?;
    let htsfile: *mut htslib::htsFile = rhtslib_hts_open(&alignment_path, b"r")?;
    Ok(htsfile)
}

pub fn rhtslib_hts_open(path: &[u8], mode: &[u8]) -> Result<*mut htslib::htsFile> {
    let cpath = ffi::CString::new(path)?;
    let c_str = ffi::CString::new(mode)?;
    let ret = unsafe { htslib::hts_open(cpath.as_ptr(), c_str.as_ptr()) };
    assert!(!ret.is_null(), "Unable to read alignment file");

    Ok(ret)
}

#[allow(clippy::not_unsafe_ptr_arg_deref)]
pub fn rhtslib_set_reference<P: AsRef<Path>>(htsfile: *mut htslib::htsFile, path: P) -> Result<()> {
    unsafe { bam::set_fai_filename(htsfile, path).context("Error setting reference file")? };
    Ok(())
}

#[allow(clippy::not_unsafe_ptr_arg_deref)]
pub fn rhtslib_fetch_by_str(
    idx: *mut htslib::hts_idx_t,
    header: *mut htslib::sam_hdr_t,
    region: &[u8],
) -> Result<*mut htslib::hts_itr_t> {
    let rstr = ffi::CString::new(region)?;
    let itr = unsafe { htslib::sam_itr_querys(idx, header, rstr.as_ptr()) };
    if itr.is_null() {
        bail!("Problem fetching reads from region")
    };

    Ok(itr)
}

pub fn rhtslib_read(
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    record: &mut Record,
) -> Option<Result<(), htsError>> {
    match rhtslib_itr_next(htsfile, itr, &mut record.inner as *mut htslib::bam1_t) {
        -1 => None,
        -2 => Some(Err(htsError::BamTruncatedRecord)),
        -4 => Some(Err(htsError::BamInvalidRecord)),
        _ => {
            // record.set_header(Rc::clone(&self.header)); // this does not seem to be necessary?
            Some(Ok(()))
        }
    }
}

#[allow(clippy::not_unsafe_ptr_arg_deref)]
pub fn rhtslib_itr_next(
    htsfile: *mut htslib::htsFile,
    itr: *mut htslib::hts_itr_t,
    record: *mut htslib::bam1_t,
) -> i32 {
    unsafe {
        htslib::hts_itr_next(
            (*htsfile).fp.bgzf,
            itr,
            record as *mut ::std::os::raw::c_void,
            htsfile as *mut ::std::os::raw::c_void,
        )
    }
}
