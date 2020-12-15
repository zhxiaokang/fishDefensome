configfile: "config/rbh_cod_gadMor3.yaml"

# The end rule for the whole workflow: what's the final output you want from this workflow?
rule end:
    input:
        cod_zf_filter = config["OUTPUTPATH"] + "/gadMor3ID_pid_gid_gname_zf_filter.tsv"

# prepare Pfam files for defensome
rule pfam:  
    input:
        pfam_dfsm_list = config["INPUTPATH"] + "/Pfam-defensome.list",
        pfamA_file = config["INPUTPATH"] + "/Pfam-A.hmm"
    output:
        pfam_dfsm_acc = config["OUTPUTPATH"] + "/Pfam-defensome.list.acc",
        pfamA_dfsm_acc = config["OUTPUTPATH"] + "/Pfam-A.hmm.defensome.acc",
        pfam_dfsm_file = config["OUTPUTPATH"] + "/Pfam-A.hmm.defensome"
    shell:
        """
        cut -f1 -d'.' {input.pfam_dfsm_list} > {output.pfam_dfsm_acc}
        grep -f {output.pfam_dfsm_acc} {input.pfamA_file} | awk '{{print $2}}' > {output.pfamA_dfsm_acc}
        hmmfetch -f {input.pfamA_file} {output.pfamA_dfsm_acc} > {output.pfam_dfsm_file}
        """

# search those defensome related Pfam in the proteome
rule hmmsearch:
    input:
        pfam_dfsm_file = config["OUTPUTPATH"] + "/Pfam-A.hmm.defensome",
        peptide = config["PEP_SEQ"]
    output:
        hmm_out = config["OUTPUTPATH"] + "/hmmsearch.tblout",
        hit_prot_id = config["OUTPUTPATH"] + "/hmmsearch.tblout.pid.uniq",
        hit_prot_seq = config["OUTPUTPATH"] + "/hmmsearch.hit.pep.fasta"
    params:
        prot_index = config["PEP_SEQ"] + ".index"
    shell:
        """
        hmmsearch -E {config[E_HMMSEARCH]} --tblout {output.hmm_out} {input.pfam_dfsm_file} {input.peptide}
        awk '/^#/ {{next}} {{print $1}}' {output.hmm_out} | sort -u > {output.hit_prot_id}
        if test -f {params.prot_index}
        then
            echo "Protein index file already exists"
        else
            fasta-make-index {input.peptide}
        fi
        fasta-fetch {input.peptide} -f {output.hit_prot_id} > {output.hit_prot_seq}
        """

# BLAST between the defensome related prot seqs and zebrafish proteome
rule blast:
    input:
        pep_zf = config["INPUTPATH"] + "/Danio_rerio.GRCz11.pep.all.fa",
        hit_prot_seq = config["OUTPUTPATH"] + "/hmmsearch.hit.pep.fasta"
    output:
        zf_as_db = config["OUTPUTPATH"] + "/blast_zf_as_db.txt",
        zf_as_query = config["OUTPUTPATH"] + "/blast_zf_as_query.txt"
    params:
        zf_db_index = config["INPUTPATH"] + "/Danio_rerio.GRCz11.pep.all.fa.pin",
        e_value = config["EVALUE_BLAST"],
        threads = config["NCORE"]
    shell:
        """
        echo "blast using zf as db"
        if test -f {params.zf_db_index}
        then
            echo "Zebrafish DB already exists"
        else
            makeblastdb -in {input.pep_zf} -dbtype prot
        fi
        blastp -db {input.pep_zf} -query {input.hit_prot_seq} -out {output.zf_as_db} -outfmt "6 qseqid sseqid evalue length pident bitscore ppos" -evalue {params.e_value} -num_threads {params.threads}
        echo "blast using zf as query"
        makeblastdb -in {input.hit_prot_seq} -dbtype prot
        blastp -db {input.hit_prot_seq} -query {input.pep_zf} -out {output.zf_as_query} -outfmt "6 qseqid sseqid evalue length pident bitscore ppos" -evalue {params.e_value} -num_threads {params.threads}
        """

# get the reciprocal best hits of the two blast results above       
rule getRBH:
    input:
        zf_as_db = config["OUTPUTPATH"] + "/blast_zf_as_db.txt",
        zf_as_query = config["OUTPUTPATH"] + "/blast_zf_as_query.txt"
    output:
        f2z = config["OUTPUTPATH"] + "/best_hit_fish_to_zf.txt",
        z2f = config["OUTPUTPATH"] + "/best_hit_zf_to_fish.txt",
        list_f2z = config["OUTPUTPATH"] + "/list_best_hit_fish_to_zf.txt",
        list_z2f = config["OUTPUTPATH"] + "/list_best_hit_zf_to_fish.txt",
        rbh = config["OUTPUTPATH"] + "/list_orthologs_fish_zf.txt"
    shell:
        """
        sort -u {input.zf_as_db} | sort -k1,1 -k3,3g -rk6,6 | awk '{{if (a[$1]=="") {{a[$1]=$3" "$6; print $0; next}}; if (a[$1]==$3" "$6) print $0}}' > {output.f2z}
        sort -u {input.zf_as_query} | sort -k1,1 -k3,3g -rk6,6 | awk '{{if (a[$1]=="") {{a[$1]=$3" "$6; print $0; next}}; if (a[$1]==$3" "$6) print $0}}' > {output.z2f}
        cut -f 1,2 {output.f2z} | sort -u > {output.list_f2z}
        cut -f 1,2 {output.z2f} | awk '{{print $2"\t"$1}}' | sort -u > {output.list_z2f}
        comm -12 {output.list_f2z} {output.list_z2f} | sort -u | awk '{{if (a[$1]=="") {{print $0; a[$1]=$2;}} else next}}' > {output.rbh}
        """

# Use the one-way best hit for those missing targets (can't find a rbh)
rule oneWay:
    input:
        f2z = config["OUTPUTPATH"] + "/best_hit_fish_to_zf.txt",
        rbh = config["OUTPUTPATH"] + "/list_orthologs_fish_zf.txt",
        hit_prot_id = config["OUTPUTPATH"] + "/hmmsearch.tblout.pid.uniq"
    output:
        rbh_fish_uniq = config["OUTPUTPATH"] + "/rbh_fish.uniq",
        miss = config["OUTPUTPATH"] + "/missing_targets.txt",
        best_hit_miss = config["OUTPUTPATH"] + "/best_hit_missing_target.txt"
    shell:
        """
        cut -f1 {input.rbh} | sort -u > {output.rbh_fish_uniq}
        comm -23 {input.hit_prot_id} {output.rbh_fish_uniq} > {output.miss}
        awk 'NR==FNR{{a[$1]=$1}} NR>FNR{{if ($1 in a) {{print $1"\t"$2}}}}' {output.miss} {input.f2z} | sort -u | awk '{{if (a[$1]=="") {{print $0; a[$1]=$2;}} else next}}' > {output.best_hit_miss}
        """

# combine the RBH and one-way best hits
rule combine:
    input:
        rbh = config["OUTPUTPATH"] + "/list_orthologs_fish_zf.txt",
        best_hit_miss = config["OUTPUTPATH"] + "/best_hit_missing_target.txt"
    output:
        bests = config["OUTPUTPATH"] + "/bests_fish_zf.txt"
    shell:
        """
        cat {input.rbh} {input.best_hit_miss} | awk -F'[\t]' '{{split($2, zf, "."); print $1"\t"zf[1]}}' > {output.bests}
        """

# retrieve gene ID and gene name of protein id from ENSEMBL, only for Zebrafish
rule retrieve:
    input:
        bests = config["OUTPUTPATH"] + "/bests_fish_zf.txt"
    output:
        zf = config["OUTPUTPATH"] + "/pid_gid_gname_zf.tsv",
        cod_zf = config["OUTPUTPATH"] + "/gadMor3ID_pid_gid_gname_zf.tsv"
    shell:
        """
        Rscript scripts/pid2gid_gname_trinity.R {input} {output.zf}
        awk 'ARGIND==1 {{zf[$1]=$2"\t"tolower($3)}} ARGIND==2 {{print $1"\t"$2"\t"zf[$2]}}' {output.zf} {input} > {output.cod_zf}
        """

# pick the zebrafish genes whose names start with the names in Marta's defensome gene list (only look at the last column!!! not the 3rd col which is the gene name of cod)
rule filter:
    input:
        fish_zf_table = config["OUTPUTPATH"] + "/gadMor3ID_pid_gid_gname_zf.tsv",
        dfsm_gene = config["INPUTPATH"] + "/defensome_genes.txt",
        dfsm_gene_cate = config["INPUTPATH"] + "/defensome_genes_category.txt",
        dfsm_letter = config["INPUTPATH"] + "/dfsm_letter.txt",
        dfsm_number = config["INPUTPATH"] + "/dfsm_number.txt",
        dfsm_any = config["INPUTPATH"] + "/dfsm_any.txt"
    output:
        cod_zf_filter = config["OUTPUTPATH"] + "/gadMor3ID_pid_gid_gname_zf_filter.tsv"
    shell:
        """
        cat {input.dfsm_gene} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($NF == LINE) print $0"\t"a[LINE]}}' {input.dfsm_gene_cate} {input.fish_zf_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        cat {input.dfsm_letter} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($NF ~ "^"LINE"[a-z].*") print $0"\t"a[LINE]}}' {input.dfsm_gene_cate} {input.fish_zf_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        cat {input.dfsm_number} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($NF ~ "^"LINE"[0-9].*") print $0"\t"a[LINE]}}' {input.dfsm_gene_cate} {input.fish_zf_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        cat {input.dfsm_any} | while read LINE
        do
            awk -v LINE="$LINE" -F'\t' 'NR==FNR{{if ($1==LINE) a[LINE]=$0}} NR>FNR{{if ($NF ~ "^"LINE"[0-9a-z].*") print $0"\t"a[LINE]}}' {input.dfsm_gene_cate} {input.fish_zf_table} >> {config[OUTPUTPATH]}/output.filter.temp
        done

        mv {config[OUTPUTPATH]}/output.filter.temp {output}
        """

