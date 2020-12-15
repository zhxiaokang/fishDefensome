# combine and merge the outputs from search_name, rbh

configfile: "config/merge_gadMor3.yaml"

rule end:
    input:
        table_merge = config["OUTPUT"]

# merge all defensome related gene names together from the 2 files
rule mergeGene:
    input:
        cds = config["CDS"],
        rbh = config["RBH"]
    output:
        dfsm_gname_cate = config["OUTPUTPATH"] + "/dfsm_gname_cate.tsv"
    shell:
        """
        awk 'NR==FNR{{print $2"\t"$4}} NR>FNR{{if(NF==6){{print $4"\t"$6}} else{{print "ERROR!! Not enough fields in file!"; exit 1}}}}' {input.cds} {input.rbh} | sort -u -k2,2 -k1,1 > {output}
        """

# combine the tables together
rule combineTable:
    input:
        cds = config["CDS"],
        rbh = config["RBH"],
        dfsm_gname_cate = config["OUTPUTPATH"] + "/dfsm_gname_cate.tsv"
    output:
        table = config["OUTPUTPATH"] + "/combine_table.tsv"
    shell:
        """
        awk -f scripts/combine_gadMor3.awk {input.cds} {input.rbh} {input.dfsm_gname_cate} > {output}
        """

# merge column cds and rbh together
rule mergeColumns:
    input:
        table = config["OUTPUTPATH"] + "/combine_table.tsv"
    output:
        table_merge = config["OUTPUT"]
    shell:
        "python scripts/merge_other_fish.py {input} {output}"

