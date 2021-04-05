# PacMonSTR-merge
Merge output from PacMonSTR

## Quick start
```
git clone https://github.com/oscarlr-TRs/PacMonSTR-merge
cd PacMonSTR-merge
python merge_pacmontr_output.py trs.bed trs.fofn merged.bed
```

## Manual
```
usage: merge_pacmontr_output.py [-h] [--stat STAT] [--min_score MIN_SCORE]
                                [--allele_diff ALLELE_DIFF]
                                trbed fofn outbed

positional arguments:
  trbed                 BED file with TRs used in PacMonSTR
  fofn                  file with paths to PacMonSTR output
  outbed                Output BED file

optional arguments:
  -h, --help            show this help message and exit
  --stat STAT           Stat to use: either avg or max
  --min_score MIN_SCORE
                        Minimum score to use
  --allele_diff ALLELE_DIFF
                        Unique allele difference
```

## Input files
### trbed
This file must be the exact bed file used in PacMonSTR with the following columns:
```
1. chrom
2. start
3. end
4. motif size
5. motif sequence
6. copies of motif in the reference
```
### fofn
This file contains the output paths from PacMonSTR. For example:
```
HG00438 /genotype_tr/HG00438/Y/4_trs.bed
HG00621 /genotype_tr/HG00621/Y/4_trs.bed
HG00673 /genotype_tr/HG00673/Y/4_trs.bed
HG00735 /genotype_tr/HG00735/Y/4_trs.bed
```
## Input optional params
### stat
Since there could be multiple reads overlapping a TR, there could be multiple copies reported by PacMonSTR. PacMonSTR-merge will either merge the max or average of those copies.

### min_score
The minimum alignment score of the expected and observed expansion sequence, and the flanks.

### allele_diff
The TR copies between sample can differ due to sequence error. For example sample A might have 3 copies and sample B might have 3.2 copies. However, they both likely have 3 copies. Therefore in order to determine the number of alleles, the alleles between the samples must differ by the `allele_diff`. If `allele_diff` is set to 1, then the alleles will differ by 1. For examples, if the alleles across the samples are `[1.5,1.8,2.2,2.5,5,10]` then the alleles would be `[1.5,2.5,5,10]`.

## Output
The output is a BED file with 13 columns:
```
1. chrom
2. start
3. end
4. motif
5. copies_in_ref
6. number_of_samples
7. number_of_genotyped_haps
8. num_alleles
9. max_alleles
10. alleles_mean
11. alleles_std
12. unique_alleles
13. all_alleles
```
