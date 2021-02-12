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
            start = line[1]
            end = line[2]
            motif_seq = line[4]
            if "/" in motif_seq:
                motif_seq = motif_seq.split("/")[0]
            copies_in_ref = line[5]
            trs[(chrom,start,end,motif_seq,copies_in_ref)] = []
    return trs

def score_filtered_copies(pacmonstr_tr,hap_copy,hap_score,min_score):
    copies = []
    prefix_scores = map(float,getattr(pacmonstr_tr,hap_score[0]).split(','))
    suffix_scores = map(float,getattr(pacmonstr_tr,hap_score[1]).split(','))
    motif_scores = map(float,getattr(pacmonstr_tr,hap_score[2]).split(','))
    hap_copies =  map(float,getattr(pacmonstr_tr,hap_copy).split(','))
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

def hap_copy(pacmonstr_tr,hap_copies,hap_scores,min_score,stat):
    copy = None
    if min_score > 0:
        all_copies = score_filtered_copies(pacmonstr_tr,hap_copies,hap_scores,min_score)
    else:
        all_copies = map(float,getattr(pacmonstr_tr,hap_copies).split(','))
    if len(all_copies) != 0:
        copy = filter_by_stat(all_copies,stat)
    return copy
        
def sample_tr_copies(pacmonstr_tr,min_score,stat):
    copies = []
    hap_copies = ["hap0_copies","hap1_copies","hap2_copies"]
    hap0_scores = ["hap0_prefix_scores","hap0_suffix_scores","hap0_motif_scores"]
    hap1_scores = ["hap1_prefix_scores","hap1_suffix_scores","hap1_motif_scores"]
    hap2_scores = ["hap2_prefix_scores","hap2_suffix_scores","hap2_motif_scores"]
    if getattr(pacmonstr_tr,"hap1_copies") == "None" and getattr(pacmonstr_tr,"hap2_copies") == "None":        
        if getattr(pacmonstr_tr,"hap0_copies") != "None":
            hap0_copy = hap_copy(pacmonstr_tr,"hap0_copies",hap0_scores,min_score,stat)
            if hap0_copy != None:
                copies = [hap0_copy]
    else:
        if getattr(pacmonstr_tr,"hap1_copies") != "None":
            hap1_copy = hap_copy(pacmonstr_tr,"hap1_copies",hap1_scores,min_score,stat)
            if hap1_copy != None:
                copies.append(hap1_copy)
        if getattr(pacmonstr_tr,"hap2_copies") != "None":
            hap2_copy = hap_copy(pacmonstr_tr,"hap2_copies",hap2_scores,min_score,stat)
            if hap2_copy != None:
                copies.append(hap2_copy)
    assert len(copies) < 3
    return copies
    
def load_sample_trs(trs,fn,sample_name,stat,min_score):
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
    Str = namedtuple('str',pacmonstr_header)
    with open(fn,'r') as fh:
        for i,line in enumerate(fh):
            line = line.rstrip().split('\t')
            if i == 0:
                if line != pacmonstr_header:
                    sys.exit("%s does not have the correct header" % fn)
                continue
            if len(line) != len(pacmonstr_header):
                print "%s is messed up" % fn
                continue
            tr = Str._make(line)
            copies = sample_tr_copies(tr,min_score,stat)
            trs[(tr.chrom,tr.start,tr.end,tr.motif,tr.copies_in_ref)].append((sample_name,copies))
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
    sample_names = []
    for sample_name,tr in trs:
        sample_added = False
        for allele in tr:
            sample_added = True
            alleles.append(allele)
        if sample_added:
            sample_names.append("%s_%s" % (sample_name,len(tr)))
    return (alleles,sample_names)

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
              "alleles_mean","alleles_std","unique_alleles","all_alleles","all_samples"]
    trs = load_trs(kwargs["trbed"])
    with open(kwargs["pacmonstr_fofn"],'r') as fofh:
        for line in fofh:
            line = line.rstrip().split('\t')
            sample_name = line[0]
            fn = line[1]
            trs = load_sample_trs(trs,fn,sample_name,kwargs["stat"],kwargs["min_score"])
    with open(kwargs["outbed"],'w') as fh:
        fh.write("%s\n" % "\t".join(header))
        for tr in trs:
            tr_copies_with_samples = trs[tr]
            alleles,sample_names = all_alleles(tr_copies_with_samples)
            tr_copies = [i[1] for i in tr_copies_with_samples]
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
                uniq_alleles = ",".join(map(str,uniq_alleles))
                alleles = ",".join(map(str,alleles))
                samples = ",".join(sample_names)
                tr_stat = [num_samples,num_haps,num_alleles,max_alleles,
                           alleles_mean,alleles_std,uniq_alleles,alleles,samples]
            else:
                tr_stat = ["."]*9
            output = list(tr) + tr_stat
            fh.write("%s\n" % "\t".join(map(str,output)))
            
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
    parser.add_argument('--allele_diff',type=int,default=1,
                        help='Unique allele difference')
    args = parser.parse_args()
    merge_trs(**vars(args))

if __name__ == '__main__':
    main()
