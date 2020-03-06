####################### RNA Seq mapping and QC pipeline for Vijay and Ay lab - Single/Paired end #######################

# Version: 1.0 (02/04/2020)

# Niu Du (ndu [at] lji.org)
# La Jolla Institute for Immunology (LJI)
# La Jolla, CA USA

# User environments
# python 3.X 

# Updates (03/06/20)
# All apps are now defined in json file
# added single-ended mode

import json
import os
import pandas as pd
from functools import reduce
from glob import glob
from QC_rule import RNA_QC
pd.set_option('display.max_colwidth', -1)

####################### Functions for reading info from pipeline generated reports #######################

def Hawk_smash(List):
    '''Flattern lists within a list '''
    return [item for sublist in List for item in sublist]

def read_fastp(sample):
    with open(f'fastp_output/{sample}_fastp.json') as json_file:
        data = json.load(json_file)
        total_reads = int(data['summary']['before_filtering']['total_reads']/2)
        filtered_reads = int(data['summary']['after_filtering']['total_reads']/2)
        filtered_reads_perc = round(100*filtered_reads/total_reads,2)
        adaptor_trimm_perc = round(100*data['adapter_cutting']['adapter_trimmed_reads']/data['summary']['before_filtering']['total_reads'],2)
        dup_rate = round(100*data['duplication']['rate'],2)
        link = f'<a href="./fastp_report/{sample}_fastp.html">seq_QC</a>'
    return [total_reads,filtered_reads,filtered_reads_perc,adaptor_trimm_perc,dup_rate,link]

def read_STAR(sample):
    with open(f'bam_aligned/{sample}/{sample}_Log.final.out') as STAR_file:
        for l in STAR_file.readlines():
            if 'Uniquely mapped reads number' in l:
                unique_reads = int(l.split('|')[-1][1:-1])
            if 'Uniquely mapped reads %' in l:
                unique_reads_perc = float(l.split('|')[-1][1:-2])    
            if 'Number of splices: Total' in l:
                spliced_reads = int(l.split('|')[-1][1:-1])
            if 'Number of splices: Annotated (sjdb)' in l:
                anno_spliced_reads = int(l.split('|')[-1][1:-1])
            if 'Number of reads unmapped: too short' in l:
                too_short_reads = int(l.split('|')[-1][1:-1])
            if '% of reads unmapped: too short' in l:
                too_short_reads_perc = float(l.split('|')[-1][1:-2])
    link = f'<a href="./bam_aligned/{sample}/{sample}_Log.final.out">map_report</a>'
    return [unique_reads,unique_reads_perc,spliced_reads,anno_spliced_reads,too_short_reads,too_short_reads_perc,link]

def read_Qualimap(sample):
    with open(f'qualimap/{sample}/rnaseq_qc_results.txt') as Qualimap_file:
        for l in Qualimap_file.readlines():
            if 'exonic' in l:
                Exonic_perc = float(l.split('=')[-1].split('(')[1][:-3])
            if 'intronic' in l:
                intronic_perc = float(l.split('=')[-1].split('(')[1][:-3])
            if 'intergenic' in l:
                intergenic_perc = float(l.split('=')[-1].split('(')[1][:-3])
            if "5' bias" in l:
                bias_5_prim = l.split('=')[-1][1:-1]
            if "3' bias" in l and "5'-3' bias" not in l:
                bias_3_prim = l.split('=')[-1][1:-1]
            if "5'-3' bias" in l:
                bias_5to3_prim = l.split('=')[-1][1:-1]
                if bias_5to3_prim == "?":
                    bias_5to3_prim = 100.0
    link = f'<a href="./qualimap/{sample}/qualimapReport.html">map_QC</a>'
    return [Exonic_perc,intronic_perc,intergenic_perc,bias_5_prim,bias_3_prim,bias_5to3_prim,link]        

def read_bamqc(sample):
    with open(f'qualimap_bamqc/{sample}/genome_results.txt') as Qualimap_file:
        for l in Qualimap_file.readlines():
            if 'mean insert size' in l:
                insert_mean = float(''.join(l.split('=')[-1][:-2].split(',')))
            if 'median insert size' in l:
                insert_median = float(''.join(l.split('=')[-1][:-1].split(',')))
    link = f'<a href="./qualimap_bamqc/{sample}/qualimapReport.html">bam_QC</a>'
    return [insert_mean,insert_median,link]    
   
def read_bamcoverage(sample): # Not used
    df = pd.read_csv(f'qualimap/{sample}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt',sep = '\t')
    return df[df.columns[1]]*100/sum(df[df.columns[1]])

def read_bw(sample,dirin,genome_version):
    relative_dir = "/".join(dirin.split("/")[3:])
    link = f'<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db={genome_version}&position=chr1&hgct_customText=track%20type=bigWig%20name={sample}%20bigDataUrl=https://informaticsdata.liai.org/NGS_analyses/ad_hoc/{relative_dir}/bed_wiggle/{sample}.bw">bigwig</a>'
    return [link]

# Here will call all previous functions and for each sample generate a list that will be inserted into the dataframe
def read_sample(sample,dirin,genome_version):
    return Hawk_smash([read_fastp(sample),read_STAR(sample),read_Qualimap(sample),read_bamqc(sample),read_bw(sample,dirin,genome_version)])


####################### Generate bash scripts for submission #######################
def gen_Submission(fastq_table,dict_conf,thread = 4,seq_type = 'Paired'):
    '''Generate bash files that are ready to be submitted to the cluster'''
    dirin = dict_conf['config']['dirin']
    bed_dir = dict_conf['config']['bed_dir']
    ref_dir = dict_conf['config']['ref_dir']
    gtf_dir = dict_conf['config']['gtf_dir']
    
    fastp = dict_conf['app']['fastp']
    STAR = dict_conf['app']['STAR']
    samtools = dict_conf['app']['samtools']
    bamCoverage = dict_conf['app']['bamCoverage']
    qualimap = dict_conf['app']['qualimap']
    
    annotation_file = dict_conf['config']['annotation_file']
    
    for sample in fastq_table.index:
        fastq_f = fastq_table.loc[sample]['fastq_f']
        fastq_r = fastq_table.loc[sample]['fastq_r']
        with open(f'Submissions/{sample}.sh','w') as f:
            # Prep
            f.write(f'#!/bin/bash\n#PBS -N {sample}\n#PBS -o {dirin}/temp/out_{sample}\n#PBS -e {dirin}/temp/err_{sample}\n#PBS -q default\n#PBS -l nodes=1:ppn={thread}\n#PBS -l mem=40gb\n#PBS -l walltime=10:00:00\ncd {dirin}/\n')
            f.write('mkdir -p fastp_output\nmkdir -p fastp_report\nmkdir -p Fastq_filtered\nmkdir -p Input\nmkdir -p counts\n')
            f.write(f'cp {ref_dir}/../{annotation_file} Input/ \n')
            # fastp
            if seq_type == 'Paired':
                f.write(f'{fastp} -w {thread} -i {fastq_f} -I {fastq_r} -o {dirin}/Fastq_filtered/{sample}_R1.fastq.gz -O {dirin}/Fastq_filtered/{sample}_R2.fastq.gz -j {dirin}/fastp_output/{sample}_fastp.json -h {dirin}/fastp_report/{sample}_fastp.html\n\n')
            elif seq_type == 'Single':
                f.write(f'{fastp} -w {thread} -i {fastq_f}  -o {dirin}/Fastq_filtered/{sample}_R1.fastq.gz -j {dirin}/fastp_output/{sample}_fastp.json -h {dirin}/fastp_report/{sample}_fastp.html\n\n')
            # STAR_MAPPING
            f.write(f'mkdir -p {dirin}/bam_aligned/{sample}\n')
            if seq_type == 'Paired':
                f.write(f'{STAR} --runThreadN {thread} --genomeDir {ref_dir} --sjdbGTFfile {gtf_dir} --readFilesIn {dirin}/Fastq_filtered/{sample}_R1.fastq.gz {dirin}/Fastq_filtered/{sample}_R2.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix {dirin}/bam_aligned/{sample}/{sample}_\n\n')
            elif seq_type == 'Single':
                f.write(f'{STAR} --runThreadN {thread} --genomeDir {ref_dir} --sjdbGTFfile {gtf_dir} --readFilesIn {dirin}/Fastq_filtered/{sample}_R1.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix {dirin}/bam_aligned/{sample}/{sample}_\n\n')
            # Export for QC
            f.write(f'mkdir -p {dirin}/bed_wiggle\n')
            f.write(f'{samtools} index {dirin}/bam_aligned/{sample}/{sample}_Aligned.sortedByCoord.out.bam\n')
            f.write(f'{bamCoverage} -p {thread}  -b {dirin}/bam_aligned/{sample}/{sample}_Aligned.sortedByCoord.out.bam -o {dirin}/bed_wiggle/{sample}.bw\n')
            f.write(f'mkdir -p {dirin}/qualimap/{sample}\n')
            if seq_type == 'Paired':
                f.write(f'{qualimap} rnaseq -pe --sorted --java-mem-size=30G -bam {dirin}/bam_aligned/{sample}/{sample}_Aligned.sortedByCoord.out.bam -gtf {gtf_dir} -outdir {dirin}/qualimap/{sample}\n')
            elif seq_type == 'Single':
                f.write(f'{qualimap} rnaseq --sorted --java-mem-size=30G -bam {dirin}/bam_aligned/{sample}/{sample}_Aligned.sortedByCoord.out.bam -gtf {gtf_dir} -outdir {dirin}/qualimap/{sample}\n')
            f.write(f'mkdir -p {dirin}/qualimap_bamqc/{sample}\n')
            f.write(f'{qualimap} bamqc -nt {thread} --java-mem-size=30G --skip-duplicated -bam {dirin}/bam_aligned/{sample}/{sample}_Aligned.sortedByCoord.out.bam -outdir {dirin}/qualimap_bamqc/{sample}\n')
    with open('2.Submissions.sh','w') as f:    
        for submission in glob('Submissions/*.sh'):        
            f.write(f'qsub {submission}\n')

    # Login herman server, go to project folder, and bash Submission.sh      

####################### Calculate TPM and export count tables, make report #######################

def merge_counts(sample_list,df_anno):
    count_tables = [pd.read_csv(f'bam_aligned/{sample}/{sample}_ReadsPerGene.out.tab',sep = '\t',skiprows = 4,header = None).set_index(0)[[1]] for sample in sorted(sample_list)]
    df_counts = reduce(lambda left,right: pd.merge(left,right,left_index = True, right_index = True), count_tables)
    df_counts.columns = sorted(sample_list)
    filter_gene_list = df_anno[df_anno["gene_type"].str.contains('tRNA|rRNA')==False].index
    df_counts_filter = df_counts.loc[filter_gene_list]
    counts = df_counts_filter.values/df_anno['Length'].loc[df_counts_filter.index].values[:,None]
    df_TPM = pd.DataFrame((1e6*counts)/(counts.sum(axis = 0)), index = filter_gene_list, columns=df_counts.columns)
    df_TPM.to_csv('counts/TPM_counts.csv')
    df_counts.to_csv('counts/raw_counts_unfiltered.csv')
    df_counts_filter.to_csv('counts/raw_counts.csv')
    return [df_counts,df_counts_filter]
    

def make_report(dict_conf,sample_list):
    dirin = dict_conf['config']['dirin']
    genome_version = dict_conf['config']['genome_version']
    annotation_file = dict_conf['config']['annotation_file']
    df_anno = pd.read_csv(f'Input/{annotation_file}').set_index('geneid')
    
    # calculate/export count table and  TPM
    [df_counts,df_counts_filter] = merge_counts(sample_list,df_anno)
    
    # Make report files
    df_QC_report = pd.DataFrame([read_sample(sample,dirin,genome_version) for sample in sample_list])
    df_QC_report.index = sample_list
    df_QC_report.columns = ['total_reads','filtered_reads','filtered_reads_perc','adaptor_trimm_perc','dup_rate','seq_QC','uniquely_mapped_reads','uniquely_mapped_reads_perc','spliced_reads','anno_spliced_reads','too_short_reads','too_short_reads_perc','map_report','exonic_perc','intronic_perc','intergenic_perc','bias_5_prim','bias_3_prim','bias_5to3_prim','map_QC','insert_mean','insert_median','bam_QC','bigwig']

    df_QC_report['STAR_counts'] =  df_counts[df_QC_report.index].sum(axis = 0)
    df_QC_report['STAR_counts_perc'] =  round(100*df_QC_report['STAR_counts']/df_QC_report['uniquely_mapped_reads'],2)
    df_QC_report['final_STAR_counts'] = df_counts_filter[df_QC_report.index].sum(axis = 0)
    df_QC_report['t_rRNA_counts'] = df_counts[df_QC_report.index].sum(axis = 0) - df_counts_filter[df_QC_report.index].sum(axis = 0)
    df_QC_report['t_rRNA_counts_perc'] = round(100*df_QC_report['t_rRNA_counts']/df_QC_report['STAR_counts'],2)

    df_genetype = df_anno[['gene_type_4']].merge(df_counts_filter,left_index = True, right_index = True).groupby('gene_type_4').sum().T
    df_QC_report['protein_coding_perc'] = (df_genetype['protein_coding']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
    df_QC_report['pseudogene_perc'] = (df_genetype['pseudogene']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
    df_QC_report['long-noncoding_perc'] = (df_genetype['long-noncoding']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
    df_QC_report['short-noncoding_perc'] = (df_genetype['short-noncoding']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
    df_QC_report['Total_genes'] = (df_counts>10).sum(axis = 0)
    df_QC_report['recommendation'] = [RNA_QC(row,dict_conf['QC_threshold']) for row in df_QC_report.itertuples()]
    
    
    df_QC_report = df_QC_report[['total_reads','filtered_reads','filtered_reads_perc','adaptor_trimm_perc','dup_rate','uniquely_mapped_reads','uniquely_mapped_reads_perc','spliced_reads','anno_spliced_reads','too_short_reads','too_short_reads_perc','exonic_perc','intronic_perc','intergenic_perc','bias_5_prim','bias_3_prim','bias_5to3_prim','STAR_counts','STAR_counts_perc','t_rRNA_counts','t_rRNA_counts_perc','protein_coding_perc','pseudogene_perc','long-noncoding_perc','short-noncoding_perc','final_STAR_counts','insert_mean','insert_median','Total_genes','seq_QC','map_report','map_QC','bam_QC','bigwig','recommendation']]

    df_QC_report = df_QC_report.sort_index()
    df_QC_report[['total_reads','filtered_reads','filtered_reads_perc','adaptor_trimm_perc','dup_rate','uniquely_mapped_reads','uniquely_mapped_reads_perc','spliced_reads','anno_spliced_reads','too_short_reads','too_short_reads_perc','exonic_perc','intronic_perc','intergenic_perc','bias_5_prim','bias_3_prim','bias_5to3_prim','STAR_counts','STAR_counts_perc','t_rRNA_counts','t_rRNA_counts_perc','protein_coding_perc','pseudogene_perc','long-noncoding_perc','short-noncoding_perc','final_STAR_counts','insert_mean','insert_median','Total_genes','recommendation']].to_csv('QC_report.csv')

    df_html = df_QC_report[['seq_QC','map_report','map_QC','bam_QC','bigwig','recommendation']]
    df_html.index.name = f'<a href="./check_QC_PCA.html">PCA_plot</a>\t\t<a href="./QC_report.csv">Download_table</a>\t\t<a href="./QC_plots">QC_plots</a>'
    df_html.to_html(f'QC_report.html',escape=False,notebook = True)
    
    df_note = pd.read_csv('/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/QC_notes.csv')
    with open(f'QC_report.html','r') as f:
        text = f.read()
        text = text.replace('</style>\n','</style>\n' + df_note.to_html(header = None,index = None).replace('border="1"','border="0"'))
    with open(f'QC_report.html','w') as f:
        f.write(text)
        
        