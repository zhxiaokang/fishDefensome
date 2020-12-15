# combine and merge the outputs from search_name, rbh

fishes = ["cod", "medaka", "stickleback", "killifish"]
fish = input("Which fish are you analysing? (Select from: cod, medaka, stickleback, killifish)\n")

if fish not in fishes:
    raise Exception("Please only select from the following fish: cod, medaka, stickleback, killifish")

if fish == "cod":
    configfile: "config/merge_cod_no_trinity.yaml"
elif fish == "medaka":
    configfile: "config/merge_mdk.yaml"
elif fish == "stickleback":
    configfile: "config/merge_stkb.yaml"
elif fish == "killifish":
    configfile: "config/merge_fdls.yaml"

rule end:
    input:
        table_merge = config["OUTPUT"]

# merge all defensome related gene names together from the 2 files
rule mergeGene:
    input:
        cdna = config["CDNA"],
        rbh = config["RBH"]
    output:
        dfsm_gname_cate = config["OUTPUTPATH"] + "/dfsm_gname_cate.tsv"
    shell:
        """
        awk 'ARGIND==1{{print $2"\t"$4}} ARGIND==2{{if(NF==8){{print $6"\t"$8}} else if(NF==7){{print $5"\t"$7}} else{{print "ERROR!! Not enough fields in file!"; exit 1}}}}' {input.cdna} {input.rbh} | sort -u -k2,2 -k1,1 > {output}
        """

# in rbh, there are some cases where fish gene names are different from zebrafish gene names
rule naughtyID:
    input:
        rbh = config["RBH"]
    output:
        naughtyID = config["OUTPUTPATH"] + "/naughtyID.tsv"
    shell:
        """
        awk 'NF==8 && $3 != $6 {{print $0}}' {input} > {output}
        """

# combine the tables together
rule combineTable:
    input:
        cdna = config["CDNA"],
        rbh = config["RBH"],
        naughtyID = config["OUTPUTPATH"] + "/naughtyID.tsv",
        dfsm_gname_cate = config["OUTPUTPATH"] + "/dfsm_gname_cate.tsv"
    output:
        table = config["OUTPUTPATH"] + "/combine_table.tsv"
    shell:
        """
        awk -f scripts/combine_other_fish.awk {input.cdna} {input.rbh} {input.naughtyID} {input.dfsm_gname_cate} > {output}
        """

# merge column cdna and rbh together
rule mergeColumns:
    input:
        table = config["OUTPUTPATH"] + "/combine_table.tsv"
    output:
        table_merge = config["OUTPUT"]
    shell:
        "python scripts/merge_other_fish.py {input} {output}"

