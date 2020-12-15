# do dea for BaP zf

Rscript do_dea.R ../input/BaP_zf/GSE64198_counts_gene_3.3_ensembl.txt ../output/gene_tmm_norm_BaP_zf_3.3_control.tsv ../output/gene_tmm_norm_BaP_zf_3.3_treat.tsv ../output/dea_BaP_zf_3.3.tsv

Rscript do_dea.R ../input/BaP_zf/GSE64198_counts_gene_96_ensembl.txt ../output/gene_tmm_norm_BaP_zf_96_control.tsv ../output/gene_tmm_norm_BaP_zf_96_treat.tsv ../output/dea_BaP_zf_96.tsv

# merge dfsm with dea
awk -F'\t' 'NR==FNR{split($1, gs, " "); g=tolower(gs[1]); a[g]=$2"\t"$5} NR>FNR{print $0"\t"a[$2]}' ../output/dea_BaP_zf_3.3.tsv ../input/defensome_gene_list_zf.txt > ../output/merge_dfsm_dea_BaP_zf_3.3.tsv

awk -F'\t' 'NR==FNR{split($1, gs, " "); g=tolower(gs[1]); a[g]=$2"\t"$5} NR>FNR{print $0"\t"a[$2]}' ../output/dea_BaP_zf_96.tsv ../input/defensome_gene_list_zf.txt > ../output/merge_dfsm_dea_BaP_zf_96.tsv