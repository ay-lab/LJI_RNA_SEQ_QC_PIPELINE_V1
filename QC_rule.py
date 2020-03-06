####################### QC rule for RNA Seq quality evaluation - Paired end #######################

# Version: N/A (02/04/2020)

# Niu Du (ndu [at] lji.org)
# La Jolla Institute for Immunology (LJI)
# La Jolla, CA USA

# User environments
# python 3.X 
# Load in threshold file from 'RNA_SEQ_QC_threshold.json'
# Please use QC parameter table at /mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/Paired_end_QC_notes.csv for reference

import json

def RNA_QC(sample_tuple,dict_threshold):
    if sample_tuple.final_STAR_counts > dict_threshold['final_STAR_counts']:
        A = True
    else:
        A = False 
    if sample_tuple.uniquely_mapped_reads_perc > dict_threshold['uniquely_mapped_reads_perc']:
        B = True
    else:
        B = False
    if sample_tuple.exonic_perc > dict_threshold['exonic_perc']:
        C = True
    else:
        C = False
    if sample_tuple.too_short_reads_perc < dict_threshold['too_short_reads_perc']:
        D = True
    else:
        D = False
    if sample_tuple.t_rRNA_counts_perc < dict_threshold['t_rRNA_counts_perc']:
        E = True
    else:
        E = False
    if sample_tuple.Total_genes  > dict_threshold['Total_genes']:
        F = True
    else:
        F = False
    if float(sample_tuple.bias_5to3_prim) < dict_threshold['bias_5to3_prim']:
        G = True
    else:
        G = False
       
    
    
    if A and B and C and D and E and F and G:
        return '1.Good'
    elif B and C and D and E and F and G:
        return '2.Reseq'
    else:
        return '3.Manual QC'