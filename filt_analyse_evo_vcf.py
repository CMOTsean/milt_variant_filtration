import pandas as pd
import filter_vcf as vcf
import numpy as np
import os


def allele_depth_ratio(s):
    s = s[0]+"/"+s[2:]
    global adr_count
    if s.startswith("0/0"):
        return s
    else:
        if vcf.read_depth_check(s, 0.25, 0.75):
            return s
        else:
            s = s.split(":")
            s[0] = "./."
            no_call = ":".join(s)
            adr_count += 1
            return no_call

def deleterious_filter(row, homo=True, stop_loss=True):
    exclusion_list = ["0/0", "./."]
    disruption_list = ["STOP-LOSS", "STOP-GAIN", "FRAMESHIFT DELETION", "FRAMESHIFT INSERTION"]
    if "|" in row.INFO:
        if row.INFO.split("|")[4] == "CDS" and row.INFO.split("|")[5] in disruption_list:
            if sum([sample[0:3] not in exclusion_list for sample in row[9:]]) > 0:
                return True
            else:
                return False
        elif row.INFO.split("|")[4] == "CDS" and row.INFO.split("|")[5] == "NONSYNONYMOUS" and row.INFO.split("|")[12] == "DELETERIOUS":
            if sum([sample[0:3] not in exclusion_list for sample in row[9:]]) > 0:
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def deleterious_count(filtered_df, strains_list):
    strain_dict = {}
    for s in strains_list:
        strain_dict[s] = {"All_genes": set(),
                          "Homo_only": set(),
                          "Het_to_hom": set()}
    for i, j in filtered_df.iterrows():
        gene = j.INFO.split("|")[2]
        for s in strains_list:
            WT = s[:3] + "WT"
            if j[s][:3] not in ["./.", "0/0"] and j[s][:3] != j[WT][:3] and j[WT][:3] != "./.":
                strain_dict[s]["All_genes"].add(gene)
                if j[s][0] == j[s][2]:
                    strain_dict[s]["Homo_only"].add(gene)
                    if j[WT][0] != j[WT][2]:
                        strain_dict[s]["Het_to_hom"].add(gene)
    return strain_dict


def gene_dict_from_strain_dict(strain_dict, filtered_df, homo=False):
    gene_dict = {}
    for i, j in filtered_df.iterrows():
        gene = j.INFO.split("|")[2]
        gene_dict[gene] = {"strains": set(),
                           "total_strains": 0}
    for g in gene_dict:
        for s in strain_dict:
            if not homo:
                if g in strain_dict[s]["All_genes"]:
                    gene_dict[g]["strains"].add(s)
                    gene_dict[g]["total_strains"] += 1
            elif homo:
                if g in strain_dict[s]["Homo_only"]:
                    gene_dict[g]["strains"].add(s)
                    gene_dict[g]["total_strains"] += 1
    for g in gene_dict.copy():
        if gene_dict[g]["total_strains"] == 0:
            del gene_dict[g]
    return gene_dict


strains = ['247A1', '247B1', '247C1', '247D1', '247D16', '247D2', '795B1', '795B16', '247E1', '247E16']
ref = "/home/sean/Desktop/Updated_ref/cpar_v3_inc_mito.fasta"
gff = "/home/sean/Desktop/Updated_ref/full_cpar_v3.gff3"
filein = "/home/sean/Desktop/Evo_analysis/SIFT4G_results/evo_samples.snps.ClusFilt.GQDP.Filtered_SIFTpredictions.vcf"
fileout = "/home/sean/Desktop/Evo_analysis/evo_filtered_deleterious.vcf"

df = pd.read_csv(filein, skiprows=vcf.skip_rows_check(filein, fileout), sep="\t")
df = df.rename(columns={"#CHROM": "CHROM"})
print(df)
df = df[vcf.repetitive_check(df, ref)]
print(df)
adr_count=0
for i in df.columns[9:]:
    df[i] = df[i].apply(np.vectorize(allele_depth_ratio))

delete_df = df[df.apply(deleterious_filter, axis=1)]
print(delete_df)
strain_dict = deleterious_count(delete_df, strains)
gene_dict = gene_dict_from_strain_dict(strain_dict, delete_df)
homo_gene_dict = gene_dict_from_strain_dict(strain_dict, delete_df, homo=True)
delete_df.to_csv(fileout, mode="a", sep="\t", header=False, index=False)

os.system(f"sed -i 's/SHIFT /SHIFT_/g' {fileout}")
