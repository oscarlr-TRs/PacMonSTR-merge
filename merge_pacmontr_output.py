#!/bin/env python
import sys
import argparse
import numpy as np
from collections import namedtuple

def load_trs(trbed):
    trs = {}
    with open(trbed,'r') as fh:
        for line in fh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        motif_seq = line[4]
        if "/" in motif_seq:
            motif_seq = motif_seq.split("/")[0]
        copies_in_ref = line[5]
        trs[(chrom,start,end,motif_seq,copies_in_ref)] = []
    return trs

def score_filtered_copies(pacmonstr_tr,hap_copy,hap_score,min_score):
    copies = []
    prefix_scores = map(float,getattr(pacmonstr_tr,hap_score[0]))
    suffix_scores = map(float,getattr(pacmonstr_tr,hap_score[1]))
    motif_scores = map(float,getattr(pacmonstr_tr,hap_score[2]))
    hap_copies =  map(float,getattr(pacmonstr_tr,hap_copy))
    for p,s,m,c in zip(prefix_scores,suffix_scores,motif_scores,hap_copies):        
        if min([p,s,m]) > min_score:
            copies.append(c)
    return copies

def filter_by_stat(copies,stat):
    assert len(copies) > 0
    if stat == "avg":
        return float(sum(copies))/len(copies)
    elif stat == "max":
        return max(copies)
    else:
        sys.exit("Choose avg or max")
    
def sample_tr_copies(pacmonstr_tr,min_score,stat):
    copies = []
    hap_copies = ["hap0_copies","hap1_copies","hap2_copies"]
    hap0_scores = [("hap0_prefix_motif_scores","hap0_suffix_motif_scores","hap0_motif_scores")]
    hap1_scores = [("hap1_prefix_motif_scores","hap1_suffix_motif_scores","hap1_motif_scores")]
    hap2_scores = [("hap2_prefix_motif_scores","hap2_suffix_motif_scores","hap2_motif_scores")]
    if getattr(pacmonstr_tr,"hap1_copies") == "." and getattr(pacmonstr_tr,"hap2_copies") == ".":        
        if getattr(pacmonstr_tr,"hap0_copies") != ".":
            copies = getattr(pacmonstr_tr,"hap0_copies")
            if min_score > 0:
                copies = score_filtered_copies(pacmonstr_tr,"hap0_copies",hap0_scores,min_score)
            if len(copies) > 0:
                copies = filter_by_stat(copies,stat)
    else:
        if getattr(pacmonstr_tr,"hap1_copies") != ".":
            copies = getattr(pacmonstr_tr,"hap1_copies")
            if min_score > 0:
                copies = score_filtered_copies(pacmonstr_tr,"hap1_copies",hap1_scores,min_score)
            if len(copies) > 0:
                copies = filter_by_stat(copies,stat)
        if getattr(pacmonstr_tr,"hap2_copies") != ".":
            if min_score > 0:
                hap2_copies = score_filtered_copies(pacmonstr_tr,"hap2_copies",hap2_scores,min_score)
                if len(hap2_copies) > 0:
                    copies += filter_by_stat(hap2_copies,stat)
            else:
                copies += filter_by_stat(getattr(pacmonstr_tr,"hap2_copies"),stat)
    return copies
    
def load_sample_trs(trs,fn,stat,min_score):
    pacmonstr_header = [
        "chrom",
        "start",
        "end",
        "motif",
        "copies_in_ref",
        "hap0_avg",
        "hap0_max",
        "hap0_copies",
        "hap0_prefix_scores",
        "hap0_suffix_scores",
        "hap0_motif_scores",
        "hap1_avg",
        "hap1_max",
        "hap1_copies",
        "hap1_prefix_scores",
        "hap1_suffix_scores",
        "hap1_motif_scores",
        "hap2_avg",
        "hap2_max",
        "hap2_copies",
        "hap2_prefix_scores",
        "hap2_suffix_scores",
        "hap2_motif_scores"
        ]
    Str = namedtupled('str',pacmonstr_header)
    with open(fn,'r') as fh:
        for i,line in enumerate(fh):
            line = line.rstrip().split('\t')
            if i == 0:
                if line != pacmonstr_header:
                    sys.exit("%s does not have the correct header" % fn)
                continue
            tr = Str._make(line)
            copies = sample_tr_copies(tr,min_score,stat)
            trs[(tr.chrom,tr.start,tr.end,tr.motif_seq,tr.copies_in_ref)].append(copies)
    return trs

def number_of_samples(trs):
    count = 0
    for tr in trs:
        if len(tr) != 0:
            count += 1
    return count

def number_of_haps(trs):
    count = 0
    for tr in trs:
        count += len(tr)
    return count

def all_alleles(trs):
    alleles = []
    for tr in trs:
        for alllele in tr:
            alleles.append(allele)
    return alleles

def unique_alleles(alleles,diff):
    smallest_allele = min(alleles)
    uniq_alleles = [smallest_allele]
    alleles = sorted(alleles)
    for allele in alleles:
        if allele > smallest_allele + diff:
            uniq_alleles.append(allele)
            smallest_allele = allele
    return uniq_alleles
        
def merge_trs(**kwargs):
    header = ["chrom","start","end","motif","copies_in_ref",
              "number_of_samples","number_of_genotyped_haps","num_alleles","max_alleles",
              "alleles_mean","alleles_std","unique_alleles","all_alleles"]
    trs = load_trs(kwargs["trbed"])
    with open(kwargs["pacmonstr_fofn"],'r') as fofh:
        for fn in fofh:
            fn = fn.rstrip()
            trs = load_sample_trs(trs,fn)
    with open(kwargs["outbed"],'w') as fh:        
        for tr in trs:
            tr_copies = trs[tr]
            alleles = all_alleles(tr_copies)
            if len(alleles) > 0:
                uniq_alleles = unique_alleles(alleles,kwargs["allele_diff"])
                
                num_samples = number_of_samples(tr_copies)
                num_haps = number_of_haps(tr_copies)
                num_alleles = len(uniq_alleles)
                max_alleles = max(uniq_alleles)
                alleles_mean = round(np.mean(np.array(alleles)),2)
                if len(alleles) > 1:
                    alleles_std = round(np.std(np.array(alleles)),2)
                else:
                    alleles_std = "."
                uniq_alleles = ",".join(uniq_alleles)
                alleles = ",".join(alleles)
                tr_stat = [num_samples,num_haps,num_alleles,max_alleles,
                           alleles_mean,alleles_std,uniq_alleles,alleles]
            else:
                tr_stat = ["."]*8
            output = list(tr) + tr_stat
            fh.write("\t".join(output))
            
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('trbed',metavar='trbed',
                        help='BED file with TRs used in PacMonSTR')
    parser.add_argument('pacmonstr_fofn',metavar='fofn',
                        help='file with paths to PacMonSTR output')
    parser.add_argument('outbed',metavar='outbed',
                        help='Output BED file')
    parser.add_argument('--stat',default="max",
                        help='Stat to use: either avg or max')
    parser.add_argument('--min_score',type=float,default=0,
                        help='Minimum score to use')
    parser.add_argument('--allele_dif',type=int,default=1,
                        help='Unique allele difference')
    args = parser.parse_args()
    merge_trs(**vars(args))

if __name__ == '__main__':
    main()
