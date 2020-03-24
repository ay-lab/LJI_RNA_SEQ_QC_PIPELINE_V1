# **LJI RNA SEQ QC PIPELINE**
* Niu Du (ndu [at] lji.org)
* La Jolla Institute for Immunology (LJI)
* La Jolla, CA USA
* Version: 1.2 (03/24/2020)

For LJI internal use, you do not need to make any changes. 

**For external users, please edit the following scripts and change the python environment to your own settings.**

Pipeline_functions/RNA_SEQ_QC_GO.py
Pipeline_functions/RNA_SEQ_sub_GO.py

change 

    #!/mnt/BioHome/ndu/anaconda3/bin/python
to 

    #!/YOUR_PYTHON_PATH 

or

    #!/usr/bin/env python

Also, please install python libraries that are listed in the requirements:

    pip install -r requirements.txt

For details please see [Document](https://ndu-ucsd.github.io/RNA_SEQ_PIPELINE/).


## Release Note

## Version: 1.2 (03/23/2020)
The following functions have been added to the pipeline:
1. Saturation curve - for detecting minimal total genes counts and minimal STAR counts of a population.
2. Outlier finder- for detecting outliers using spearman correlation and minimal threshold of 0.6.
3. Note - specific reason why a perticular sample failed the QC test is now available in the Note column.

The following rule has been implemented:
1. Good sample if : Star counts > 3x10^6 reads (paired) or minimal threshold detected, whoever is higher; uniquely mapped reads > 60%, too short reads < 50%; exotic reads > 50%; STAR counts in total uniquely mapped reads > 80%; t/rRNA < 10%;Total genes per sample > 5000 or minimal threshold detected, whoever is higher;  0.95 < 5’-3’ bias < 1.25; median insert size > 150 bp.
2. Reseq sample if: Only Star counts failed while other criteria were good.
3. Manual QC: Any other criterion did not met minimal requirement.

## Version: 1.1 (02/25/2020)

The following functions have been added to the pipeline:
1. PCA plot - for QC and quick visualization of metadata (if provided); sample name is visible when hovering over scatter point. To change category please use the dropdown list on top left corner. 
2. USCS genome browser interface - additional bigwig link has been added to the report webpage that allows automatic upload of bw file for mapping visualization. (Thanks to Cristian for contribution to the function.) 

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
4. To modify tools (ver only), ~~the command lines are located in STEP 2 of RNA_SEQ_Func.py~~,make changes in the json file.
5. Current setting is optimized based on discussion with Seq team; in case of change to RNA_SEQ_Func.py, the pipeline need to be debugged thoroughly
6. Add editing log below in this document

### Note 
In case of deep sequencing depth (e.g. more than 50 million reads) you might experience problem of insurficient memory for the QC step; If that happened, do not rush to the 3rd step, stop here, go to the Submissions folder and find corresponding sh file, and make necessary changes to the sh file. You do not need to change the pipeline. 
