import pandas as pd
import numpy as np
import filter_vcf as vcf
#filein = "/home/sean/Desktop/all_samples.snps.ClusFilt.GQDP.Filtered.vcf"
#fileout = "/home/sean/Desktop/all_samples.CustomFilt.vcf"
ref = "/home/sean/Desktop/Updated_ref/cpar_v3_inc_mito.fasta"

##Test
filein = "/home/sean/Desktop/genotyped_samples.Excluded.vcf"
fileout = "/home/sean/Desktop/genotyped_samples.CustomFilt.vcf"


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

adr_count=0
df = pd.read_csv(filein, skiprows=vcf.skip_rows_check(filein, fileout), sep="\t")
df = df.rename(columns={"#CHROM": "CHROM"})

print(df.shape)
rep_check = vcf.repetitive_check(df, ref)
df = df.iloc[rep_check]
print(df.shape)
#test = df.iloc[4883:4885,:]
for i in df.columns[9:]:
    df[i] = df[i].apply(np.vectorize(allele_depth_ratio))
print(adr_count)

no_call_truth_list=[]
for index, r in df.iterrows():
    if sum(r.str.startswith('./.', na=False)) == (df.shape[1]-9):
        no_call_truth_list.append(False)
    else:
        no_call_truth_list.append(True)
df = df.iloc[no_call_truth_list]
print(df.shape)

###COMMENT THIS OUT
index_list=[]
for i, r in df.iterrows():
    if (r["247"][0:3] != r["247WT"][0:3]) and (r["247"][0:3] != "./." and r["247WT"][0:3] != "./."):
        index_list.append(i)
for i, r in df.iterrows():
    if (r["795"][0:3] != r["795WT"][0:3]) and (r["795"][0:3] != "./." and r["795WT"][0:3] != "./."):
        index_list.append(i)


df = df.loc[index_list]
df.to_csv(fileout, mode="a", sep="\t", header=False, index=False)






