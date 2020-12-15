
import pandas as pd
import sys
import copy

file_input = sys.argv[1]
file_output = sys.argv[2]

df = pd.read_csv(file_input, sep = "\t", header = None)

cdna_id = list(df.iloc[:, 2])
cdna_id = [x for x in cdna_id if str(x) != "nan"]

dic_cate = {}
dic_origin = {}
for row in range(df.shape[0]):
    cate = df.iloc[row, 0]
    key = df.iloc[row, 1]
    cdna = df.iloc[row, 2]
    rbh = df.iloc[row, 3]
    if str(key) != "nan":
        gene = key
        dic_cate[gene] = cate
        dic_origin[gene] = (list(), list())
        if str(cdna) != "nan":
            dic_origin[gene][0].append(cdna)
        if str(rbh) != "nan":
            dic_origin[gene][1].append(rbh)
    else:
        if str(cdna) != "nan":
            dic_origin[gene][0].append(cdna)
        if str(rbh) != "nan":
            dic_origin[gene][1].append(rbh)

dic_merge = {}  # in the merged dic, only one list is needed
for gene in dic_origin:
    dic_merge[gene] = copy.deepcopy(dic_origin[gene][0])  # every gene ID in cdna is correct, so always keep them
    for id_rbh in dic_origin[gene][1]:
        if id_rbh not in cdna_id:
            dic_merge[gene].append(id_rbh)

df_merge = pd.DataFrame(columns=["category", "gene_name", "cdna", "rbh", "merge"])
index_df_merge = 0
for gene in dic_cate:
    cdna = dic_origin[gene][0]
    rbh = dic_origin[gene][1]
    merge = dic_merge[gene]
    len_cdna = len(cdna)
    len_rbh = len(rbh)
    len_merge = len(merge)
    index_gene = 0
    while index_gene < max(len_cdna, len_rbh, len_merge):
        if index_gene == 0:
            cell_cate = dic_cate[gene]
            cell_gene = gene
        else:
            cell_cate = ""
            cell_gene = ""
            
        cell_cdna = cdna[index_gene] if index_gene < len_cdna else ""
        cell_rbh = rbh[index_gene] if index_gene < len_rbh else ""
        cell_merge = merge[index_gene] if index_gene < len_merge else ""
        
        df_merge.loc[index_df_merge] = [cell_cate, cell_gene, cell_cdna, cell_rbh, cell_merge]
        
        index_gene += 1
        index_df_merge += 1

df_merge.to_csv(file_output, index=False, sep="\t")

