#!/usr/bin/env python3
import argparse
import json
import os
import time

from cyvcf2 import VCF
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

VERSION = "1.0.0"

SUPPORTED = (
    ".json",
    ".eps", 
    ".jpeg",
    ".jpg", 
    ".pdf",
    ".pgf",
    ".png",
    ".ps",
    ".raw",
    ".rgba",
    ".svg",
    ".svgz",
    ".tif",
    ".tiff",
    ".webp",
)

VCF_SKIP_TAGS = ["CNMISSING", "CNZERO", "DPZERO"]

DESCRIPTION="\
description:\n\
    Generate an overview of normalised depth of coverage values from ConSTRain VCF\n\
    output. Normalised depth of coverage is determined by dividing the number of\n\
    reads mapping to a repeat locus by the copy number (i.e., the number of alleles\n\
    that exists for this locus in the genome). This information is useful to rerun\n\
    ConSTRain on a VCF file with new --min-norm-depth and --max-norm-depth values\n\
    based on the observed depth distribution to filter out outlier loci.\
" 
def parse_cla():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        epilog=f"Script version v{VERSION}",
    )


    parser.add_argument(
        "-v", "--vcf", type=str, required=True,
        help="VCF file output by ConSTRain for which to generate sequencing depth distribution plot"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="File handle where output plot will be saved. Filetype of output will be based on the file extension. \
            If extension is json, a json file with the bounds will be generated. Otherwise, a histogram will be generated."
    )
    parser.add_argument(
        "-a", "--alpha", type=float, required=True,
        help="Used to calculate bounds for loci to include. Lower bound will be at quantile alpha/2, upper bound at quantile 1 - (alpha/2). \
            (NOTE: lower bound will always be at least 1., as it is impossible to estimate a genotype with fewer reads than there are alleles)"
    )
    parser.add_argument(
        "--include_mononuc", action="store_true",
        help="Include mononucleotide repeats when determining bounds. \
            (NOTE: mononucleotide repeats will still be included in the final \
                output plot no matter what, they just won't be considered for determining the bounds \
                    when --include_mononuc is not set)"
    )

    return parser.parse_args()

def df_from_vcf(vcf_file: str) -> pd.DataFrame:
    vcf = VCF(vcf_file)
    if len(vcf.samples) != 1:
        raise RuntimeError("this script currently only supports analysing VCF files with exactly one sample")
    df = {
        "period": [],
        "copy_number": [],
        "depth": [],
    }

    for variant in vcf:
        if any([variant.format("FT")[0] == tag for tag in VCF_SKIP_TAGS]):
            continue
        df["period"].append(variant.INFO.get("PERIOD"))
        df["copy_number"].append(variant.format("CN")[0][0])
        df["depth"].append(variant.format("DP")[0][0])
    
    df = pd.DataFrame(df).assign(depth_norm = lambda x: x["depth"] / x["copy_number"])
    return df

def main():    
    args = parse_cla()
    ext = os.path.splitext(args.output)[-1]
    if ext not in SUPPORTED:
        raise ValueError(f"File extension {ext} is not a supported output format. Use one of {SUPPORTED}")

    print(f"Parsing VCF file {args.vcf}")
    start = time.time()    
    df = df_from_vcf(args.vcf).dropna()
    print(f"Read VCF file in {time.time() - start:.2f} seconds")    

    if args.include_mononuc:
        lower = max(1., df["depth_norm"].quantile(q=args.alpha/2, interpolation="nearest"))
        upper = df["depth_norm"].quantile(q=1 - (args.alpha/2), interpolation="nearest")
    else:        
        lower = max(1., df.query("period > 1")["depth_norm"].quantile(q=args.alpha/2, interpolation="nearest"))
        upper = df.query("period > 1")["depth_norm"].quantile(q=1 - (args.alpha/2), interpolation="nearest")
    n_within = df.query(f"depth_norm >= {lower} and depth_norm <= {upper}").shape[0]
    
    if ext == ".json":
        print("Writing to json")
        output = {
            "lower": lower,
            "upper": upper,
            "n_within": n_within,
            "n": df.shape[0],
        }
        with open(args.output, 'w') as f:
            json.dump(output, f, indent=4)
    else:
        print("Generating histogram")
        start = time.time()

        xmin = max(0, lower - 15)
        xmax = upper + 15
        df_plot = df.query(f"depth_norm >= {xmin} and depth_norm <= {xmax}")

        sns.set_context("poster")
        fig = plt.figure(dpi=300)    
        ax = sns.histplot(
            df_plot,
            x = "depth_norm",
            discrete=True,
            stat="proportion",
            color="grey",
        )    

        _ = ax.set(
            xlabel="Depth / CN",
            xlim=(xmin, xmax),
            ylabel="Proportion of STR loci",
            title=f"Lower bound: {lower}, upper bound: {upper}\nLoci in range: {n_within}/{df.shape[0]} ({n_within/df.shape[0]*100:.2f}%)",
        )

        _ = ax.vlines(
            x=[lower, upper],
            ymin=0, 
            ymax=ax.get_ylim()[1], 
            color="black",
        )
        print(f"Generated histogram in {time.time() - start:.2f} seconds")

        print(f"Saving plot to {args.output}")
        plt.savefig(args.output, bbox_inches="tight")
        plt.close()

if __name__ == "__main__":
    main()
