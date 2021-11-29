#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:41:15 2020

@author: sean
"""

import pandas as pd
import numpy as np

def sequence_check(chrom, query_pos, infile, rep_list):
    for i in infile:
        if i.startswith(chrom):
            sequence = i.replace("\n", "")[len(chrom)+query_pos:len(chrom)+query_pos+8]
            rev_sequence = i.replace("\n", "")[len(chrom)+query_pos-9:len(chrom)+query_pos-1]
    if sequence not in rep_list and rev_sequence not in rep_list:
        return True
    else:
        return False
        
    
def repetitive_check(df, ref_file):
    with open(ref_file, "r") as ref_file:
        file = ref_file.read().split(">")
    rep_list=["TTTTTTTT",
              "AAAAAAAA",
              "CCCCCCCC",
              "GGGGGGGG",
              "TATATATA",
              "TCTCTCTC",
              "TGTGTGTG",
              "GAGAGAGA",
              "GCGCGCGC",
              "GTGTGTGT",
              "CACACACA",
              "CTCTCTCT",
              "CGCGCGCG",
              "ATATATAT",
              "AGAGAGAG",
              "ACACACAC"]
    truth_list=[]
    for row in df.itertuples():
        truth_list.append(sequence_check(row.CHROM, row.POS, file, rep_list))
    return truth_list

def load_bad_genes(filename):
    with open(filename, "r") as bad_in:
        bad_genes_list = [l.rstrip() for l in bad_in.readlines() if not l.startswith("#")]
    return bad_genes_list

def skip_rows_check(filename, fileout):
    count=0
    with open(filename, "r") as vcf_in, open(fileout, "w") as vcf_out:
        for i in vcf_in.readlines():
            if i.startswith("#"):
                count += 1
                vcf_out.write(i)
            else:
                break
    return count-1

def load_gff_to_df(filename):
    commented_lines=[]    
    colnames=["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    with open(filename, "r") as gff_in:
        for k, l in enumerate(gff_in):
            if l.startswith("#"):
                commented_lines.append(k)
    gff_df = pd.read_csv(filename, skiprows=commented_lines, sep="\t", index_col=False, names=colnames, usecols=colnames)
    gff_df["start"] = pd.to_numeric(gff_df["start"])
    gff_df["end"] = pd.to_numeric(gff_df["end"])
    return gff_df

def check_position(chromo, query_pos, gff_df, name=False):
    query_pos = int(query_pos)
    gene_name = ""
    for row in gff_df[gff_df["feature"]=="CDS"].itertuples():
        if chromo == row.sequence:
            if query_pos in range(row.start, row.end+1):
                gene_name = row.attributes.split("Name=",1)[1].split(";")[0]
                break
    if name==True:
        if gene_name=="":
            return None
        else:
            return gene_name
    else:
        if gene_name=="":
            return False
        else:
            return True

def read_depth_check(info, het_lim, hom_lim):
    gt, ad, dp = info.split(":")[0:3]
    if dp == "." or dp == "0":
        return False
    if gt == "1/1":
        hom = True
    elif gt == "./.":
        return False
    else:
        hom = False
    #ad_ref, ad_alt = map(int, ad.split(","))
    ad_ref = int(ad.split(",")[0])
    ad_alt = int(ad.split(",")[1])
    dp = int(dp)
    if hom == False:
        if (ad_ref/dp) > het_lim and (ad_alt/dp) > het_lim:
            return True
        else:
            return False
    else:
        if (ad_alt/dp) > hom_lim:
            return True
        else:
            return False
        
        
###Main execution        
strain="CLIB214"
mode = "indel"
if mode == "indel":
    vcf = "/home/sean/Desktop/Fang_data/VCFs/"+strain+".Selected.indels.vcf"
else:
    vcf = "/home/sean/Desktop/Fang_data/VCFs/"+strain+".Selected.snps.vcf"
    
def main():
    gff_df = load_gff_to_df("/home/sean/LabDrive/MSKCC-RTA3/Parapsilosis/Annotation/full_cpar.gff3")
    bad_genes_list = load_bad_genes("/home/sean/LabDrive/MSKCC-RTA3/Parapsilosis/Annotation/bad_genes.txt")
    ref_genome = "/home/sean/LabDrive/MSKCC-RTA3/Parapsilosis/Annotation/cpar.fasta"
    df= pd.read_csv(vcf, skiprows=skip_rows_check(vcf), sep="\t")
    #df= pd.read_csv("/home/sean/Desktop/Fang_data/VCFs/62-65_indels_wo_CLIB.vcf", skiprows=42, sep="\t")
    df = df.rename(columns={"#CHROM": "CHROM", strain: "STATS"})

    count=0
    index_list=[]



    ###For Manual Checking
    #if mode == "snp":
    #    for row in df.itertuples():
    #        gene_name = check_position(row.CHROM, row.POS, gff_df, name=True)
    #        if read_depth_check(row.STATS) and gene_name:
    #            if gene_name not in bad_genes_list:
    #                #print(gene_name)
    #                count += 1
    #                index_list.append(row.Index)
    #
    #elif mode == "indel":
    #    for row in df.itertuples():
    #        gene_name = check_position(row.CHROM, row.POS, gff_df, name=True)
    #        if gene_name:
    #            if gene_name not in bad_genes_list:
    #                count += 1
    #                index_list.append(row.Index)

    ###For LOH
    if mode == "snp":
        for row in df.itertuples():
            gene_name = check_position(row.CHROM, row.POS, gff_df, name=True)
            if read_depth_check(row.STATS):
                if gene_name not in bad_genes_list:
                    #print(gene_name)
                    count += 1
                    index_list.append(row.Index)

    elif mode == "indel":
        for row in df.itertuples():
            gene_name = check_position(row.CHROM, row.POS, gff_df, name=True)
            if gene_name not in bad_genes_list:
                count += 1
                index_list.append(row.Index)

    filtered_df = df.iloc[index_list]
    rep_check = repetitive_check(filtered_df, "/home/sean/LabDrive/MSKCC-RTA3/Parapsilosis/Annotation/cpar.fasta")
    filtered_df = filtered_df.iloc[rep_check]
    #filtered_df.to_csv("/home/sean/Desktop/"+strain+"_FullFilt_"+mode+".vcf", mode="a", sep="\t", header=False, index=False)
    filtered_df.to_csv("/home/sean/Desktop/"+strain+"_LOHFilt_"+mode+".vcf", mode="a", sep="\t", header=False, index=False)

#def main():   
if __name__ == "__main__":
    main()