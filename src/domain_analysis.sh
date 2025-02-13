#!/bin/bash

dir_proj=/mnt/chaelab/rachelle/zzOtherzz/XiaoMei/HDstress
dir_domain=${dir_proj}/results/domain
dir_hmm=/mnt/chaelab/rachelle/data/hmm_profiles/profiles
dir_brara=/mnt/chaelab/rachelle/data/Brapa

#######################
##  NB-ARC DOMAIN
#######################

dir_nbarc=${dir_domain}/NB-ARC
mkdir -p ${dir_nbarc}

grp=Col0-Brapa3-Brapa3.5-R500-Z1-Cai2021

f_tmp=${dir_nbarc}/tmp.txt
domains=( TIR RX-CC_like RPW8 NB-ARC )
for domain in ${domains[@]}; do
    aln=/mnt/chaelab/rachelle/nlr_survey/results/aln/apt_TCRN/${domain}.${grp}.pep.mafft.fasta
    grep '>' ${aln} | sed 's/>//' > ${f_tmp}
    python3 -c "import re, sys; sys.path.append('/mnt/chaelab/rachelle/src'); import data_manip; seqids = data_manip.splitlines('${f_tmp}'); tids = [re.search('[^|]+(?=\|${domain}\||\|revcomp\|${domain}\|)', seqid).group(0) for seqid in seqids]; fout = open('${dir_nbarc}/${domain}.${grp}.id.txt', 'w+'); fout.write('\n'.join(['\t'.join([tids[i], seqids[i]]) for i in range(len(seqids))])); fout.close()"
done
rm ${f_tmp}


######################
##  PAN_2 DOMAIN
######################

########### find PAN_2 ################
hmm_pan=${dir_hmm}/interproscan-5.56-89.0_PAN2/interproscan-5.56-89.0_PAN2.hmm
dir_pan=${dir_domain}/PAN_2
dir_pan_hmm=${dir_pan}/hmmer
mkdir -p ${dir_pan} ${dir_pan_hmm}

## find in v3.5
fa_brara35=${dir_brara}/v3.5/Brara_Chiifu_V3.5_pep.fa
hmmsearch --tblout ${dir_pan_hmm}/Brapa_v3.5.PAN_2.hmmer.out \
          -A ${dir_pan_hmm}/Brapa_v3.5.PAN_2.hmmer.sto \
          ${hmm_pan} ${fa_brara35} > /dev/null

## find in Z1
fa_z1=${dir_brara}/Z1/GCA_900412535.3_brapa_z1_v2_protein.faa
hmmsearch --tblout ${dir_pan_hmm}/Brapa_Z1.PAN_2.hmmer.out \
          -A ${dir_pan_hmm}/Brapa_Z1.PAN_2.hmmer.sto \
          ${hmm_pan} ${fa_z1} > /dev/null

# ## find in R500 (eh. this one doesn't have pep fa :()
# fa_r500=${dir_brara}/Z1/GCA_900412535.3_brapa_z1_v2_protein.faa
# phmmer ${fa} ${hmm_pan} -o ${dir_pan_hmm}/Brapa_Z1.PAN_2.hmmer

## find in in cai2021
dir_cai=${dir_brara}/cai2021
for fa in ${dir_cai}/*.pep.fa; do
    accid=$(basename ${fa%*.pep.fa})
    out_pref=Brapa_${accid}.PAN_2.hmmer
    hmmsearch --tblout ${dir_pan_hmm}/${out_pref}.out \
              -A ${dir_pan_hmm}/${out_pref}.sto \
              ${hmm_pan} ${fa} > /dev/null
done

## find in a thaliana
fa_ath=/mnt/chaelab/shared/genomes/TAIR10/fasta/TAIR10_pep_20101214.fasta
hmmsearch --tblout ${dir_pan_hmm}/Arath_Col-0.PAN_2.hmmer.out \
          -A ${dir_pan_hmm}/Arath_Col-0.PAN_2.hmmer.sto \
          ${hmm_pan} ${fa_ath} > /dev/null

########### convert stockholm to fasta ################
dir_pan_seq=${dir_pan}/seq
mkdir -p ${dir_pan_seq}

grp=Col0-Brapa3.5-Z1-Cai2021

fa_pan=${dir_pan_seq}/${grp}.PAN_2.pep.fasta
f_tmp=${dir_pan_seq}/tmp.fasta
python3 -c "from Bio import SeqIO; import os, itertools; records = itertools.chain(*[[r for r in SeqIO.parse('${dir_pan_hmm}/' + fname, 'stockholm')] for fname in os.listdir('${dir_pan_hmm}') if os.path.splitext(fname)[1] == '.sto']); SeqIO.write(records, '${f_tmp}', 'fasta')"
/mnt/chaelab/rachelle/src/degap.py ${dir_pan_seq}/tmp.fasta ${fa_pan}
rm ${f_tmp}

########### align #############
dir_pan_aln=${dir_pan}/aln
mkdir -p ${dir_pan_aln}

grp=Col0-Brapa3.5-Z1-Cai2021

aln_pan=${dir_pan_aln}/$(basename ${fa_pan%*.fasta}).mafft.fasta
mafft ${fa_pan} > ${aln_pan}

########## tree ###############
dir_pan_tree=${dir_pan}/tree
mkdir -p ${dir_pan_tree}

grp=Col0-Brapa3.5-Z1-Cai2021

nwk_pan=${dir_pan_tree}/$(basename ${aln_pan%*.fasta}).nwk
/mnt/chaelab/rachelle/programmes/FastTree < ${aln_pan} > ${nwk_pan}

## move to desktop for plotting




######################
##  S_locus_glycop DOMAIN
######################

########### find S_locus_glycop ################
hmm_slg=${dir_hmm}/interproscan-5.56-89.0_Slocusglycop/interproscan-5.56-89.0_Slocusglycop.hmm
dir_slg=${dir_domain}/S_locus_glycop
dir_slg_hmm=${dir_slg}/hmmer
mkdir -p ${dir_slg} ${dir_slg_hmm}

## find in v3.5
fa_brara35=${dir_brara}/v3.5/Brara_Chiifu_V3.5_pep.fa
hmmsearch --tblout ${dir_slg_hmm}/Brapa_v3.5.S_locus_glycop.hmmer.out \
          -A ${dir_slg_hmm}/Brapa_v3.5.S_locus_glycop.hmmer.sto \
          ${hmm_slg} ${fa_brara35} > /dev/null

## find in Z1
fa_z1=${dir_brara}/Z1/GCA_900412535.3_brapa_z1_v2_protein.faa
hmmsearch --tblout ${dir_slg_hmm}/Brapa_Z1.S_locus_glycop.hmmer.out \
          -A ${dir_slg_hmm}/Brapa_Z1.S_locus_glycop.hmmer.sto \
          ${hmm_slg} ${fa_z1} > /dev/null

# ## find in R500 (eh. this one doesn't have pep fa :()
# fa_r500=${dir_brara}/Z1/GCA_900412535.3_brapa_z1_v2_protein.faa
# phmmer ${fa} ${hmm_slg} -o ${dir_slg_hmm}/Brapa_Z1.S_locus_glycop.hmmer

## find in in cai2021
dir_cai=${dir_brara}/cai2021
for fa in ${dir_cai}/*.pep.fa; do
    accid=$(basename ${fa%*.pep.fa})
    out_pref=Brapa_${accid}.S_locus_glycop.hmmer
    hmmsearch --tblout ${dir_slg_hmm}/${out_pref}.out \
              -A ${dir_slg_hmm}/${out_pref}.sto \
              ${hmm_slg} ${fa} > /dev/null
done

## find in a thaliana
fa_ath=/mnt/chaelab/shared/genomes/TAIR10/fasta/TAIR10_pep_20101214.fasta
hmmsearch --tblout ${dir_slg_hmm}/Arath_Col-0.S_locus_glycop.hmmer.out \
          -A ${dir_slg_hmm}/Arath_Col-0.S_locus_glycop.hmmer.sto \
          ${hmm_slg} ${fa_ath} > /dev/null

########### convert stockholm to fasta ################
dir_slg_seq=${dir_slg}/seq
mkdir -p ${dir_slg_seq}

grp=Col0-Brapa3.5-Z1-Cai2021

fa_slg=${dir_slg_seq}/${grp}.S_locus_glycop.pep.fasta
f_tmp=${dir_slg_seq}/tmp.fasta
python3 -c "from Bio import SeqIO; import os, itertools; records = itertools.chain(*[[r for r in SeqIO.parse('${dir_slg_hmm}/' + fname, 'stockholm')] for fname in os.listdir('${dir_slg_hmm}') if os.path.splitext(fname)[1] == '.sto']); SeqIO.write(records, '${f_tmp}', 'fasta')"
/mnt/chaelab/rachelle/src/degap.py ${dir_slg_seq}/tmp.fasta ${fa_slg}
rm ${f_tmp}

########### align #############
dir_slg_aln=${dir_slg}/aln
mkdir -p ${dir_slg_aln}

grp=Col0-Brapa3.5-Z1-Cai2021

aln_slg=${dir_slg_aln}/$(basename ${fa_slg%*.fasta}).mafft.fasta
mafft ${fa_slg} > ${aln_slg}

########## tree ###############
dir_slg_tree=${dir_slg}/tree
mkdir -p ${dir_slg_tree}

grp=Col0-Brapa3.5-Z1-Cai2021

nwk_slg=${dir_slg_tree}/$(basename ${aln_slg%*.fasta}).nwk
/mnt/chaelab/rachelle/programmes/FastTree < ${aln_slg} > ${nwk_slg}

## move to desktop for plotting
