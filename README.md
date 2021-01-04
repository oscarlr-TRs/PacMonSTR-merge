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
