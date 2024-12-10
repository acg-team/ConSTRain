#!/usr/bin/env python3
import argparse

from cyvcf2 import VCF
import numpy as np
import pandas as pd

VERSION = "1.0.0"

DESCRIPTION="\
description:\n\
    Create CSV file based on ConSTRain VCF output. CSV file will have six columns:\n\
        str_id:         {chromosome id}_{start position} (0-based).\n\
        copy_number:    the number of alleles that exists for this locus in the genome.\n\
        frequencies:    string representation of Python dict. Keys are allele length, values are observed frequencies.\n\
        genotype:       string representation of Python list. List the allele lengths of the inferred genotype.\n\
        depth:          the number of reads that mapped to this locus\n\
        depth_norm:     depth divided by copy_number.\
" 

def parse_cla():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=DESCRIPTION,
            epilog=f"Script version v{VERSION}",
        )

    parser.add_argument(
        "-v", "--vcf", type=str, required=True,
        help="VCF file output by ConSTRain from which to create a CSV file" 
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="File path where the CSV file should be written"
    )

    return parser.parse_args()

def df_from_vcf(vcf_file: str) -> pd.DataFrame:
    vcf = VCF(vcf_file)
    if len(vcf.samples) != 1:
        raise RuntimeError("this script currently only supports analysing VCF files with exactly one sample")
    df = {
        "str_id": [],
        "copy_number": [],        
        "frequencies": [],
        "genotype": [],
        "depth": [],
    }

    for variant in vcf:
        df["str_id"].append(f"{variant.CHROM}_{variant.POS - 1}")
        df = parse_constrain_format_field(df, variant)
    
    df = pd.DataFrame(df).assign(depth_norm = lambda x: x["depth"] / x["copy_number"])
    return df

def parse_constrain_format_field(df: dict, variant) -> dict:
    try:
        df["copy_number"].append(variant.format("CN")[0][0])
    except TypeError:
        df["copy_number"].append(np.nan)

    try:
        df["depth"].append(variant.format("DP")[0][0])
    except TypeError:
        df["depth"].append(np.nan)
    
    try:
        frequencies = variant.format("FREQS")[0]
        freq_dict = dict()
        for i in frequencies.split("|"):
            i = i.split(",")
            freq_dict[int(i[0])] = int(i[1])
        df["frequencies"].append(freq_dict)            
    except (TypeError, IndexError):
        df["frequencies"].append(np.nan)

    try:
        genotypes = variant.format("REPLEN")[0]
        if genotypes == ".":
            raise ValueError
        genotypes = [int(i) for i in genotypes.split(",")]
        df["genotype"].append(genotypes)
    except (TypeError, ValueError):
        df["genotype"].append(np.nan)
    
    return df

def main():
    args = parse_cla()

    df = df_from_vcf(args.vcf)

    df.to_csv(args.output, index=False, header=True)

if __name__ == "__main__":
    main()
