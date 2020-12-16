awk -F'\t' 'NR==FNR{a[$1]=$0} NR>FNR{print $2"\t"$1"\t"a[$3]}' ../input/expression.txt ../input/dfsm_gene_cate_id_zf.txt | sort > ../output/dfsm_tpm.txt
