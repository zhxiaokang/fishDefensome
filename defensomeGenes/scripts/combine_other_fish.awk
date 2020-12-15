# awk file to combine outputs from search_name and rbh

ARGIND==1 {name=$2; id=$1;
if (!(name in cdna))
    cdna[name]=id
else if (cdna[name]!~id)
    cdna[name]=cdna[name]"\t"id
}
ARGIND==2 {if(NF==8){name=$6; id=$2} else if(NF==7){name=$5; id=$2} else{print "ERROR!! Not enough fields in file!"; exit 1}
if (!(name in rbh))
    rbh[name]=id
else if (rbh[name]!~id)
    rbh[name]=rbh[name]"\t"id
}
ARGIND==3 {name=$6; cod_name=$3;
if (!(name in naughtyID))
    naughtyID[name]=cod_name
else if (naughtyID[name]!~cod_name)
    naughtyID[name]=naughtyID[name]"\t"cod_name
}
ARGIND==4 {name=$1; cate=$2;
len_cdna=split(cdna[name],list_cdna,"\t");
len_rbh=split(rbh[name],list_rbh,"\t");
len_naughtyID=split(naughtyID[name],list_naughtyID,"\t");
if (len_cdna < len_rbh)
    len_max = len_rbh
else
    len_max = len_cdna
if (len_max < len_naughtyID)
    len_max = len_naughtyID
print cate"\t"name"\t"list_cdna[1]"\t"list_rbh[1]"\t"list_naughtyID[1]
if (len_max > 1){
    for (i=2; i<=len_max; i++){
        print "\t\t"list_cdna[i]"\t"list_rbh[i]"\t"list_naughtyID[i]
    }
}
}
