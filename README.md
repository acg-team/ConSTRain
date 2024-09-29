# ConSTRain: copy number guided STR allele inference

ConSTRain is a short tandem repeat (STR) variant caller that can account for copy number variants (CNVs), aneuploidies, and polyploid genomes.
It is a very accurate and fast method, needing less than 20 minutes to genotype >1.7 million STRs in an 100X alignment of human whole-genome sequencing reads.
ConSTRain generates all possible allele length distributions, and returns the one that best matches the observed allele length distribution for each STR locus.

<!-- ![](images/method_overview.png) -->

<img src="./images/method_overview.png" alt="method overview" width=1000 height="auto">

## Installation
### Building from source
In order to build ConSTRain from the source code, you will need to [install the Rust programming language](https://www.rust-lang.org/tools/install).
Then, clone this repository and build ConSTRain as follows:

```bash
git clone https://github.com/acg-team/ConSTRain.git
cd ConSTRain/ConSTRain
cargo build --release --bin ConSTRain
```

The binary will be created at `./target/release/ConSTRain`.

### Possible issues
#### libclang is not found
Building `hts-sys` fails with the following error: 
```bash
Unable to find libclang: "couldn't find any valid shared libraries matching: ['libclang.so', 'libclang-*.so', 'libclang.so.*' 'libclang-*.so.*'], set the LIBCLANG_PATH environment variable to a path where one of these files can be found (invalid: [])"
```

This is caused because `Bindgen` requires `libclang`. The problem can be solved by installing `clang` (see instructions [here](https://rust-lang.github.io/rust-bindgen/requirements.html)) or, if it is already installed, providing the location of the directory containing `libclang.so` through the `LIBCLANG_PATH` environment variable:
```bash
LIBCLANG_PATH=/path/to/clang/lib cargo build --release --bin ConSTRain
```

#### System header files are not found
Building `hts-sys` fails with an error like this:
```bash
./htslib/htslib/hts.h:31:10: fatal error: 'stddef.h' file not found
```

Caused because some system header file (e.g., `stddef.h`) is not found. Address this by providing the location of the directory containing system header files through the `C_INCLUDE_PATH` environment variable:
```bash
C_INCLUDE_PATH=/path/to/include cargo build --release --bin ConSTRain
```

## Running ConSTRain
Calling STR variants from aligned sequencing reads is done via the `ConSTRain alignment` subcommand.
To run, `ConSTRain alignment` needs a BED file specifiying locations of STRs in the reference genome, a JSON file encoding the target chromosomes and their ploidies, and an alignment in SAM/BAM/CRAM format.
A basic ConSTRain invocation looks like this:

```bash
ConSTRain alignment -a <ALIGNMENT> -k <KARYOTYPE> -r <REPEATS>
```

If you know there are CNVs in the sample you're analysing, ConSTRain can account for these if you supply the affected regions as a BED file via the `--cnvs` flag:

```bash
ConSTRain alignment -a <ALIGNMENT> -k <KARYOTYPE> -r <REPEATS> --cnvs <CNVs>
```

### Re-running ConSTRain on a VCF file
Once an alignment has been analysed using `ConSTRain alignment` the generated VCF file can be used to re-run ConSTRain.
This is useful if you decide you want to use different filtering parameters, or maybe you have generated new CNV information that was not available when you first analysed a sample.
Rather than having to analyse the alignment again, you can use the `ConSTRain vcf` subcommand to re-estimate STR genotypes.
This is possible because ConSTRain includes the observed allele length distribution for each STR as a FORMAT field in the VCF file.
Running `ConSTRain vcf` is typically a matter of a few seconds.

The basic invocation of `ConSTRain vcf` is:

```bash
ConSTRain vcf -v <VCF> -k <KARYOTYPE> -r <REPEATS>
```

Again, you can add CNV information via the `--cnv` flag:

```bash
ConSTRain vcf -v <VCF> -k <KARYOTYPE> -r <REPEATS> --cnvs <CNVs>
```


## Input files
ConSTRain expects alignments to adhere to file specifications maintained by the [GA4GH](https://samtools.github.io/hts-specs/), and outputs variant calls in VCF format.
However, specific formats are expected for other input files, which are outlined here. 

*   Karyotype: the organism's karyotype should be specified as a JSON file mapping the names of chromosomes to their ploidies.
E.g., `{"chr1": 2, ... "chrX": 2, "chrY: 0"}` for a human female sample, `{"chr1": 2, ... "chrX": 1, "chrY: 1"}` for a human male.
It is critical that chromosome names in this file match chromosome names in the alignment file exactly.
*   Repeats: the location of repeats in the reference genome should be provided as a BED3+2 file, with the two extra columns indicating the repeat unit length and the repeat unit.
E.g., `chr5 21004   21014   2   AT` specifies a dinucleotide repeat with sequence `ATATATATAT` starting at position 21004 on chromosome 5.
*   Copy number variants: CNVs can be specified as a BED3+1 file, with the extra column indicating the (absolute!) copy number of the region affected by the CNV.
E.g., `chr5 106680   106750   3`.
