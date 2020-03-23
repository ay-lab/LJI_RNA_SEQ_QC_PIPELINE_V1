**LJI Sequencing core RNA-Seq workflow**
======

---
### **Prerequisites**

Install all required python libraries; we recommend using conda to manage your vitual environment if outside LJI environment.

	>conda create -n RNA_SEQ python=3.7
	>conda activate RNA_SEQ
	>which python

	>/anaconda3/envs/RNA_SEQ/bin/python

for starters at LJI you can use:

	>/mnt/BioHome/ndu/anaconda3/bin/python

---

### **RNA QC logistics**


**- Run samples**

I. Make a clone of this repository to your own directory. If you do not have required packages pre-installed, please download and install them on the cluster server following the links in the main page.

II. Edit the both .sh files to enure the link to the python scripts are correct; for example, if your cloned directory is 	

	>/home/YOURNAME/RNA_SEQ_PIPELINE

change the excution line in 1.RNA\_Seq\_Mapping.sh to 	

	>/anaconda3/envs/RNA_SEQ/bin/python /home/YOURNAME/RNA_SEQ_PIPELINE/Pipeline_functions/RNA_SEQ_sub_GO.py
	-c conf_RNA_Seq.json -i fastq_table.csv -n 4 -seq "Paired"

III. The ```RNA_SEQ_sub_GO.py```  script takes 4 inputs: ```-c``` - config file, ```-i``` fastq location table, ```-n``` cpus, ```-seq``` single ended "Single" or pair ended "Paired" input fastq type. For best performance please use the default cpu setting if possible.


* Example of fastq table; note that the names in the sample ID will be used for naming all of the down stream files.

| sample ID  | fastq_f  | fastq_r  |
|---|---|---|
|  sample 1 |  Fastq\_input/sample\_1_R1.fastq.gz | Fastq\_input/sample\_1_R2.fastq.gz  |
|  sample 2 |  Fastq\_input/sample\_2_R1.fastq.gz | Fastq\_input/sample\_2_R2.fastq.gz  |
|  sample 3 |  Fastq\_input/sample\_3_R1.fastq.gz | Fastq\_input/sample\_3_R2.fastq.gz  |



IV. Modify json file based on your settings . In most cases you will only need to edit the "dirin" and "metadata_dir" in "config", for changing reference genome please go to [Reference genome](), and for changing QC criteria please go to [QC parameters]().



#### Architecture of the json configuration file
<img src = './img/json_file.png'>

- optional: Example of metadata file; if the number of samples in metadata file is less than the total number of samples processed by the pipeline, only those in the metadata file will be shown in the PCA plot. The value of metadata_dir can be left empty (set to " "), and only QC results will be used in the PCA plot.

| sample ID  | Sex  | Disease  |
|---|---|---|
|  sample 1 |  Male | yes  |
|  sample 2 |  Female | No  |
|  sample 3 |  Male | No  |



V. Execute ```>qsub 1.RNA_Seq_Mapping.sh``` on the server; this run will create necessary folders and generate bash files in the 'Submission' folder for running the pipelines in parallel, and a '2.Submissions.sh' file in the working directory. Once finished please go to the Submission folder and double check one or two files and make sure the settings were properly configured in the scripts if you are not certain.

VI.  Execute ```>bash 2.Submissions.sh``` on the server; this step will qsub all bash files in the Submission folder. If your sample number is more than that the cluster could handle, considering separate them into different batches.

VII. Once all submissions have been completed, execute ```>qsub 3.RNA_Seq_QC.sh``` to generate QC report and plots. The 'QC_report.html' will be generated inside the working directory with links to other plots and tables.

---
**- Update the sequencing run log ([VD_Vijay&Ay_Bioinfo_log](https://docs.google.com/spreadsheets/d/1WPmbCofdCcxCRrpx6tPM78g34nTWHJTLZe7ySV3v-no/edit#gid=0))**
* note: Fill the 'Seq QC location' column with the QC report link using this format ```=HYPERLINK("https://informaticsdata.liai.org/NGS_analyses/ad_hoc/Groups/YOUR_DIR.../QC_report.html","report")```

At this point you need to contact the person in charge of the project to manually check samples in category 3 (manual QC) following the [QC tutorial](). Once finished the 'failed' samples together with the category 2 (reseq) samples will be resequenced by the sequence team; you need to notify them.  

---
**- Update libraries**

Once the re-sequencing step is finished, you will get the updated fastq files from the sequence team (or bcl files where you can extract the fastq files from). If the sample was only resequenced due to low sequencing depth without using new library, merge the old and new fastq files together and qsub the corresponding bash file in the Submission folder, and then redo step VII. If an new library was used for sequencing or the sequencing run had issues, delete the old fastq file and do the above steps.

Usually the QC step can be done within one iteration, however in case when multiple iterations are required, the same principle guidelines can be applied.
