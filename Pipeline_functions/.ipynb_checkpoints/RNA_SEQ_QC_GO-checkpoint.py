#!/mnt/BioHome/ndu/anaconda3/bin/python

__author__ = "Niu Du"
__email__ = "ndu@lji.org"

import argparse,json
from glob import glob

from RNA_SEQ_Func import merge_counts,make_report
from RNA_SEQ_QC_plots import bam_plot,QC_plot,Plot_3D,RNA_QC_spearman,soft_threshold


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--conf', type=str)
    opts = parser.parse_args()
    conf_file = opts.conf
    
    dict_conf = read_json(conf_file)
    dirin = dict_conf['config']['dirin']
    sample_list = [x.split('/')[-1].split('.')[0] for x in glob(f'{dirin}/bed_wiggle/*')]
    make_report(dict_conf,sample_list)
    
    bam_plot()
    vis_parameters = ['final_STAR_counts']
    QC_plot(vis_parameters,dict_conf['QC_threshold'])
    
    vis_parameters = ['uniquely_mapped_reads_perc','exonic_perc']
    QC_plot(vis_parameters,dict_conf['QC_threshold'])

    vis_parameters = ['too_short_reads_perc','t_rRNA_counts_perc']
    QC_plot(vis_parameters,dict_conf['QC_threshold'])
    
    vis_parameters = ['Total_genes']
    QC_plot(vis_parameters,dict_conf['QC_threshold'])
    
    vis_parameters = ['bias_5to3_prim']
    QC_plot(vis_parameters,dict_conf['QC_threshold'])
    
    vis_parameters = ['insert_median']
    QC_plot(vis_parameters,dict_conf['QC_threshold'])
    
    meta_input = dict_conf['config']['metadata_dir']
    Plot_3D(meta_input,n_comps = 10) 
    
    RNA_QC_spearman()
    
    soft_threshold(dict_conf)
    
def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf

if __name__ == "__main__":
    main()