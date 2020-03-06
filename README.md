# **LJI_RNA_SEQ_PIPELINE - AY and VIJAY lab**
* Niu Du (ndu [at] lji.org)
* La Jolla Institute for Immunology (LJI)
* La Jolla, CA USA
* Version: 1.1 (02/25/2020)

## User tutorial 
Please following the tutorial below:

1. copy files '1.RNA_Seq_Mapping.sh', '3.RNA_Seq_QC.sh', and 'conf_RNA_Seq.json' to your own working directory 
2. make changes to the conf_RNA_Seq.json file according to your settings
    * Please do not modify QC_threshold unless confirmed with sequence team and the bioinformatican who is in charge of RNA QC
3. Create your own fastq_table.csv file in your working dir follow the example in this folder; the address for fastq files can be relative (to your working dir) or absolute dir. 
4. At your working directory, do the following steps
    4.1 qsub 1.RNA_Seq_Mapping.sh; this will generate all bash files for submission to the grid, and generate a sh file for batch submission
    4.2 wait for 4.1 to finish, and bash 2.Submissions.sh; All sh files will be submitted to the grid
    4.3 wait for all tasks in 4.2 to finish, and qsub 3.RNA_Seq_QC.sh; this step will generate QC report for sequence team to review
5. Once the workflow is completed, sending report link to Seqteam, and
    5.1 in case of rerun, merge new fastq file with old ones and rerun coresponding bash file(s) in the Submission folder
    5.2 in case of redo-lib, replace orignal fastq file with new one and rerun coresponding bash file(s) in the Submission folder
    5.3 Once finished, qsub 3. RNA_Seq_QC.sh to generate new report for next round of QC.

## Version: 1.1 (02/25/2020)

The following functions have been added to the pipeline:
1. PCA plot - for QC and quick visualization of metadata (if provided); sample name is visible when hovering over scatter point. To change category please use the dropdown list on top left corner. 
Example: https://informaticsdata.liai.org/NGS_analyses/ad_hoc/Groups/vd-vijay/ndu/SeAs_TH2_Resting/check_QC_PCA.html (internal use)
2. USCS genome browser interface - additional bigwig link has been added to the report webpage that allows automatic upload of bw file for mapping visualization. (Thanks to Cristian for contribution to the function.) 
Example: https://informaticsdata.liai.org/NGS_analyses/ad_hoc/Groups/vd-vijay/ndu/SeAs_TH2_Resting/QC_report.html (internal use)

The following rule has been implemented:
1. Good sample if : Star counts > 3x10^6 reads (paired), uniquely mapped reads > 60%, too short reads < 50%, exotic reads > 50%, STAR counts in total uniquely mapped reads > 80%, t/rRNA < 10%,Total genes per sample > 5000,  5’-3’ bias < 1.25.
2. Reseq sample if: Only Star counts < 3x10^6 while other criteria were good.
3. Manual QC: Any other criterion did not met minimal requirement.

Minor bug fixed:
Fixed problem while individual user may not have proper package installed (Thanks Cristian and Vicente for testing)
Fixed problem no plot been generated while not good samples exist
Fixed problem when ? Symbol showed up in  5’-3’ bias report
Removed html_dir from Jason file  


## Version: 1.0 (02/04/2020)



####### Instruction for making changes to the pipeline #######
1. In case of major change to the pipeline, please make a copy of all files with in /mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/RNA_SEQ_Script/
2. To change the QC rule, the logic for ranking sequence is in quality QC_rule.py, and the threshold is located in the conf_RNA_Seq.json file
3. To add more plot functions, the module for each kind of plot is in RNA_SEQ_plots.py, and add excution step in RNA_SEQ_QC_GO.py
4. To modify tools (ver only), the command lines are located in STEP 2 of RNA_SEQ_Func.py
5. Current setting is optimized based on discussion with Seq team; in case of change to RNA_SEQ_Func.py, the pipeline need to be debugged thoroughly
6. Add editing log below in this document

### Note 02/11/2020
In case of deep sequencing depth (e.g. more than 50 million reads) you might encount problem of insurficient memory for the QC step; If that happened, do not rush to the 3rd step, stop here, go to the Submissions folder and find corresponding sh file, and make necessary changes to the sh file. You do not need to change the pipeline. 


### Note 02/12/2020
Fixed the error that when only 1 category show up no figs were plotted
