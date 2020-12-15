import sys

configfile: "config/cds_gadMor3.yaml"

rule end:
    input:
        name_id_table = config["OUTPUTPATH"] + "/name_id_table.txt"
#dfsm_genes_category = config["OUTPUTPATH"] + "/defensome_genes_category.txt"

rule extract:
    input:
        cdna = config["CDS"] 
    output:
        name_id_table = config["OUTPUTPATH"] + "/name_id_table.txt"
    shell:
        """
        awk '/^>/ && $2~/gene/ {{split($1, id, ">"); split($2, name, "[]=]"); print tolower(name[2])"\t"id[2]}}' {input} | sort -u > {output}
        """

rule search:
    input:
        dfsm_genes_all = config["INPUTPATH"] + "/defensome_genes.txt",
        dfsm_genes_category_all = config["INPUTPATH"] + "/defensome_genes_category.txt",
        name_id_table = config["OUTPUTPATH"] + "/name_id_table.txt",
        dfsm_letter = config["INPUTPATH"] + "/dfsm_letter.txt",
        dfsm_number = config["INPUTPATH"] + "/dfsm_number.txt",
        dfsm_any = config["INPUTPATH"] + "/dfsm_any.txt"
    output:
        dfsm_genes_category = config["OUTPUTPATH"] + "/defensome_genes_category.txt"
    shell:
        """
        cat {input.dfsm_genes_all} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($1 == LINE) print $2"\t"$1"\t"a[LINE]}}' {input.dfsm_genes_category_all} {input.name_id_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        cat {input.dfsm_letter} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($1 ~ "^"LINE"[a-z].*") print $2"\t"$1"\t"a[LINE]}}' {input.dfsm_genes_category_all} {input.name_id_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done
        
        cat {input.dfsm_number} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($1 ~ "^"LINE"[0-9].*") print $2"\t"$1"\t"a[LINE]}}' {input.dfsm_genes_category_all} {input.name_id_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        cat {input.dfsm_any} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($1 ~ "^"LINE"[a-z0-9].*") print $2"\t"$1"\t"a[LINE]}}' {input.dfsm_genes_category_all} {input.name_id_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        mv {config[OUTPUTPATH]}/output.filter.temp {output}
        """
