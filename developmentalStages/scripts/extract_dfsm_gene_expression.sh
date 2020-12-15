for file in ../input*_gene_norm.tsv
do
	awk 'NR==FNR{a[$1]=$0} NR>FNR{print $2"\t"$1"\t"a[$3]}' $file ../input/dfsm_gene_cate_id_stkb.txt | sort > ../output/dfsm_$file
done