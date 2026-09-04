#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 02.2 — Clean the MEGAN6 taxonomy export
# ══════════════════════════════════════════════════════════════════════════
# Converts the MEGAN6 export (OTU_ID + raw taxon path) into a clean
# 8-rank taxonomy path, filling "unknown/unclassified/unassigned/NA" gaps
# with the last known taxonomic rank rather than leaving them blank.
#
# INPUT:  assigned_seqs-MARES_e20-ex.txt  (MEGAN6 export, renamed)
# OUTPUT: taxonomy_edited_8ranks_e20.txt

cp -av assigned_seqs-MARES_e20-ex.txt assigned_seqs-MEGAN_e20-ex.raw.txt

awk -F'\t' 'BEGIN{
  OFS="\t";
  print "OTU_ID","TaxonPath"
}
{
  otu=$1;
  tax=$2;

  if (tax=="" || tax ~ /^[[:space:]]*$/) tax="d__unassigned;";

  gsub(/\r/,"",tax);
  gsub(/[[:space:]]+$/,"",tax);

  if (tax !~ /;$/) tax = tax ";";

  n=split(tax,a,";");
  last="";
  out="";

  for(i=1;i<=n;i++){
    if(a[i]=="" ) continue;

    if(a[i] ~ /^[dkpcofgs]__$/){
      prefix=substr(a[i],1,3);
      if(last!="") a[i]=prefix substr(last,4);
    }

    if(a[i] ~ /(unclassified|unknown|unassigned|NA|none|null)/i){
      prefix=substr(a[i],1,3);
      if(last!="") a[i]=prefix substr(last,4);
    }

    if(a[i]!="") last=a[i];
    out = out a[i] ";";
  }

  if(out=="" ) out="d__unassigned;";
  print otu,out;
}' assigned_seqs-MARES_e20-ex.txt > taxonomy_edited_8ranks_e20.txt

wc -l taxonomy_edited_8ranks_e20.txt
awk -F'\t' 'NR>1 && $2=="d__unassigned;"{c++} END{print "unassigned:",c}' taxonomy_edited_8ranks_e20.txt
