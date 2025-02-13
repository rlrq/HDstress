#!/bin/bash

dir_proj=/mnt/chaelab/rachelle/zzOtherzz/XiaoMei/HDstress

## filter orthofam data for brassica and arabidopsis genes only
f_orthofam_full=${dir_proj}/data/orthofam/genefamily_data.ORTHOFAM.csv
f_orthofam_filt=${dir_proj}/data/orthofam/genefamily_data.ORTHOFAM.ath-bra.tsv
f_orthofam_bra_id=${dir_proj}/data/orthofam/genefamily_data.ORTHOFAM.bra.id
awk '$2=="ath" || $2=="bra"' ${f_orthofam_full} > ${f_orthofam_filt}
awk '$2=="bra"' ${f_orthofam_filt} | cut -f1 | sort | uniq > ${f_orthofam_bra_id}

## map go terms to ath genes in orthofam

