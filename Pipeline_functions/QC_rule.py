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



def RNA_QC(sample_tuple,dict_threshold,count_X_low,count_Y_low):
    
    minimal_counts = dict_threshold['minimal_counts']
    
    if minimal_counts == 'fixed':# Reseq based fixed total STAR count 
        if sample_tuple.final_STAR_counts > dict_threshold['final_STAR_counts']:
            A = True
        else:
            A = False 
    elif minimal_counts == 'perc': # Reseq based on % of gene recovery rate if the threshold is higher than the fixed value
        if sample_tuple.final_STAR_counts > max(count_X_low,dict_threshold['final_STAR_counts']):
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
    # 
    if sample_tuple.Total_genes  < max(dict_threshold['Total_genes'][0],min(dict_threshold['Total_genes'][1],count_Y_low)):
        F = False
    else:
        F = True
    if (float(sample_tuple.bias_5to3_prim) < dict_threshold['bias_5to3_prim'][1])&(float(sample_tuple.bias_5to3_prim) > dict_threshold['bias_5to3_prim'][0]): # In rear cases the bias value is small; also not accepted
        G = True
    else:
        G = False
    if (sample_tuple.insert_median > dict_threshold['insert_median'][0]):# and (sample_tuple.insert_median < dict_threshold['insert_median'][1]):
        H = True
    else:
        H = False
    
    
    
    if A and B and C and D and E and F and G and H:
        return '1.Good','None'
    elif B and C and D and E and F and G and H:
        return '2.Reseq','counts='+str(int(sample_tuple.final_STAR_counts))+'threshold='+str(int(max(count_X_low,dict_threshold['final_STAR_counts'])))
    else:
        return '3.Manual QC','/'.join(name for item,name in zip([A,B,C,D,E,F,G,H],['final_STAR_counts','uniquely_mapped_reads_perc','exonic_perc','too_short_reads_perc','t_rRNA_counts_perc','Total_genes','bias_5to3_prim','insert_median']) if item == False)