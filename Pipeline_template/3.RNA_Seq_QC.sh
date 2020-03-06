#!/bin/bash

#PBS -N RNA_Seq_QC
#PBS -o RNA_Seq_QC_out
#PBS -e RNA_Seq_QC_err
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -l walltime=1:00:00

cd $PBS_O_WORKDIR
mkdir -p QC_plots
/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/RNA_SEQ_Script/RNA_SEQ_QC_GO.py -c conf_RNA_Seq.json
