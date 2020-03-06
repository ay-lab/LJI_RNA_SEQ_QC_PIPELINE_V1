#!/bin/bash

#PBS -N RNA_Seq_Mapping
#PBS -o RNA_Seq_Mapping_out
#PBS -e RNA_Seq_Mapping_err
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=1:00:00

cd $PBS_O_WORKDIR
mkdir -p Submissions
mkdir -p temp
/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/RNA_SEQ_Script/RNA_SEQ_sub_GO.py -c conf_RNA_Seq.json -i fastq_table.csv -n 4 -seq "Paired"
# -seq "Single"
