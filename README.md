# allele_stats.py: a program to calculate U20, U50, and Q95 statistics from VCF files.


## Overview 

**Goal**: To calculate various allele frequency statistics described in [(Racimo et al. 2017)](https://doi.org/10.1093/molbev/msw216). These stats include U20, U50, and Q95, which are all window-based statistics used to detect candidate sites of adaptive introgression across the genome.

**Statistic definitions**: 

- *U* statistics: "*U<sub>A,B,C</sub>*(*w,x,y*): Number of sites in which any allele is at a frequency lower than w in panel A, higher than x in panel B and equal to y in panel C" (Table 1 in Racimo et al. 2017). Currently, this script considers only w=1% and y=100% (as shown in Fig. 1 in Racimo et al. 2017). In practice, I fixed sites in panel A to the reference allele and sites in panel C to the alternate allele, ensuring that allele frequencies meet requirements of w and y in the above definition. For *U20*, x=20 and for *U50*, x=50.

- *Q95* statistic:
"*Q95<sub>A,B,C</sub>*(*w,y*): 95% quantile of the distribution of derived allele frequencies 
in panel B, for sites where the derived allele is at a frequency lower 
than w in panel A and equal to y in panel C" (Table 1 in Racimo et al. 2017). Currently, this script considers only w=1% and y=100% (as shown in Fig. 1 in Racimo et al. 2017). In practice, I fixed sites in panel A to the reference allele and sites in panel C to the alternate allele, ensuring that allele frequencies meet requirements of w and y in the above definition.

**Input**: 
- VCF file (.vcf)
- population key file (.txt)
- genomic windows file (.bed)

**Output**:
- File 1 (*allele_stats_by_site.csv*) is organized like a VCF file, but only includes sites for which population A is fixed for the reference allele and population C is fixed for the alternate allele. The script adds four additional columns that provide allele counts and frequencies used in the calculation of the U20, U50, and Q95 statistics (which are window-based). These allele counts and frequencies are pretty interesting on their own and provide specific site information (unlike the window-based statistics in the second file).

- File 2 (*allele_stats_by_window.csv*) is organized by the windows specified in the genomic window file. Each window includes the calculated value for U20, U50, Q95, and the number of informative sites in the window.

**System requirements**: 
- Currently, this scripts loads the whole VCF file into memory, so you need a lot of memory for large VCF files. During testing, I ran VCF files of up to 25 GB in size, which finished in ~20 minutes using ~155 GB of memory. If users want to reduce memory requirements, they can divide VCF files into smaller files (say, by chromosome) and process that way.
- Program developed with Python v3.11.4. Likely compatible with multiple versions, but not yet tested.

**Dependencies**:
- argparse module
- pandas module (v2.1.4)
- collections module
- numpy module (v1.26.1)

**Notes**: 

- This script filters out sites with >= 50% missingness in any of the three specified populations to help ensure data quality.
- To properly use this program it is important to identify the ancestral allele (similar logic to ABBA-BABA tests). Our best approximation of this ancestral allele is often the allele in a near outgroup. Thus, ideally, your reference genome is a species or population just outside the clade of interest. When this is the case, we can straightforwardly polarize the alleles of population A by fixing them to alternate state "1" and those of population C by fixing them to the ancestral state "0". Sites with the opposite polarization are not as informative since shared alleles between populations B and C may indicate shared ancestry more than adaptive introgression. 
- Populations A and B are sister. Population C is sister to populations A + B. Reference is outgroup of populations A + B + C.


## Setup

- Check out the example files in this repository
- Double check that your VCF file has the expected nine intro columns (from 'CHROM' through 
'FORMAT'). If you are missing the FORMAT column, add it so that the column indexing is not thrown off.
- Filter VCF to include only biallelic SNPs. Missing data is ok (see note above).
- Use command `grep -n "#CHROM" filename.vcf` to identify the line number of the VCF file 
header line. Subtract one and use that value for the `--skipRows` flag (see usage below).
- Create a tab-delimited population file (.txt) with two columns: 1) Population and 2) Individual. This file serves to map your individuals to the three populations under consideration. Values in the Population column should consist of a single word followed by an underscore and number (e.g., blueRobin_1). The word (=the population name) serves to group populations and the number keeps track of the individuals in that population. Also, use the word to identify populations when running the program (e.g., `--popA blueRobin`)
- Create a tab-delimited genomic window file (.bed) that specifies the windows over which you want to estimate U20, U50, and Q95. This file has no column headers and three columns with 1) chromosome names, 2) window start position, and 3) window end position. In order to get this genomics window file, I suggest the following steps:
    + create a genome file with two columns (no headers) containing 1) chromosome names and 2) their length in bases. The command below creates this genome file from the reference genome fasta index file (.fasta.fai): `awk -v OFS='\t' {'print $1,$2'} genome.fasta.fai > genome_file.txt`
    + Use this genome file in conjuction with bedtools `makewindows` to create the genomic window file. Here are two examples of commands to make windows (the first produces non-overlapping 10kb windows and the second produces sliding 10kb windows with 1kb overlap): \
    `bedtools makewindows -g genome_file.txt -w 10000 > windows.bed` \
    `bedtools makewindows -g genome_file.txt -w 10000 -s 9000 > windows.bed`

- Move your input files into the same directory as allele_stats.py and run the program!

## Usage

`python allele_stats.py --vcfFile filename.vcf --skipRows n
--windowFile windows.bed --popKey popKey.txt --popA popA_name
--popB popB_name --popC popCname`

## Command to run example files

`python allele_stats.py --vcfFile test.vcf --skipRows 81 --windowFile windows.bed --popKey popKey.txt --popA belem --popB xingu --popC tapajos`



## Citation

[![DOI](https://zenodo.org/badge/746977070.svg)](https://zenodo.org/doi/10.5281/zenodo.10553643)



**Moncrieff, A.E.** 2024. Allele_stats.py: a program to calculate U20, U50, and Q95 statistics from VCF files. DOI:10.5281/zenodo.10553644
