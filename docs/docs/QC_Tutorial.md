# **QC tutorial**

---
### Step 1: overview

To get an idea about the overall quality of RNASeq result, you can refer to the [recommendation and outlier](#Outlier) columns in the main QC report page. In general, you should find most of your samples labeled with 1.Good and No in those two columns, respectively. If not, there might be  critical issues in sample prep or sequencing steps that need to be diagnosed separately. The links to **PCA plot**, **QC parameters** and **QC plots** are located on the top left corner of the report table, while the sample specific QC reports are located on the right side of their corresponding sample IDs.

<img src='img/QC_report.png'>



### Step 2: Examine PCA

The interactive PCA plot allows you to quickly inspect meta information of each sample. If no meta table was provided in the json file, you can still check the recommendation and Outlier information in respect to each sample. To select meta information, click on the dropdown menu and select the target item. You can make specific category invisible by click its symbol at the bottom left corner.

<img src='img/PCA_plot.png' width = 600>



### Step 3: Examine QC table

The QC table includes all important QC parameters and should be your main source of information for manual QC. We used fairly stringent criteria for ranking the sample quality, therefore it is pretty safe to ignore a sample marked as a good sample and not a outlier (in most cases a good sample should not be a outlier). If a category 2 or 3 sample was not in the good list but close enough due to the stringent criteria, you might not need to re-sequence it. Besides these scenarios, the sample need to be resequenced or the library should be redone. If not possible then the sample should be eliminated. Below is the dictionary of QC parameters used.   

|Parameter                 |Note                                                                                                                   |
|--------------------------|-----------------------------------------------------------------------------------------------------------------------|
|total_reads               |All paired end reads; each pair is considered one read                                                                 |
|filtered_reads            |fastp filtered reads with low complexity and low quality reads removed                                                 |
|filtered_reads_perc       |percentage of remaining reads after fastq QC                                                                           |
|adaptor_trimm_perc        |percentage of reads that have adaptor content and have been trimmed                                                    |
|dup_rate                  |percentage of fastq sequence that have duplicated reads                                                                |
|uniquely_mapped_reads     |reads only mapped to one location of the genome model                                                                  |
|uniquely_mapped_reads_perc|percentage of uniquely mapped reads                                                                                    |
|spliced_reads             |total number of splicing events in each read                                                                           |
|anno_spliced_reads        |Splicing known in splice junction database                                                                             |
|too_short_reads           |the overlap between reads and genome is less than the minimal set level; in STAR the default is 60% of seqeuence length|
|too_short_reads_perc      |percentage of those too short reads                                                                                    |
|exonic_perc               |percentage of reads mapped to exonic regions                                                                           |
|intronic_perc             |percentage of reads mapped to intronic regions                                                                         |
|intergenic_perc           |percentage of reads mapped to intergenic regions                                                                       |
|bias_5_prim               |mean expression of 5' divided by mean expression of transcript                                                         |
|bias_3_prim               |mean expression of 3' divided by mean expression of transcript                                                         |
|bias_5to3_prim            |the ratio between both 5' and 3' biases                                                                                |
|STAR_counts               |STAR generated count table sum by samples                                                                              |
|STAR_counts_perc          |STAR counts in total uniquely mapped reads                                                                             |
|t_rRNA_counts_perc        |percentage of tRNA and rRNA in the total STAR counts                                                                   |
|protein_coding_perc       |percentage of protein coding gene counts in the total STAR counts                                                      |
|pseudogene_perc           |percentage of  pseudogene counts in the total STAR counts                                                              |
|long-noncoding_perc       |percentage of long-noncoding mRNA counts in the total STAR counts                                                      |
|short-noncoding_perc      |percentage of short-noncoding mRNA counts in the total STAR counts                                                     |
|final_STAR_counts         |STAR count counts sum by sample excluding tRNA and rRNA counts                                                         |
|insert_mean               | mean insert size of paired end reads                                                                                  |
|insert_median             | median insert size of paired end reads                                                                                |
|Total_genes                 | Total number of genes identified after mapping        |

### Step 4: Examine QC plots
The QC plots folder contains three types of figures that can assist your diagnoses of RNASeq sample quality in step 3.
- The first type is a simple visualization of QC parameters with dash lines mark where the set thresholds are compare to the real data (QC parameter in file names);

<img src='img/type_1_example.png' width = 1000>

In these figures, Good sample points were condensed to the left side for allowing more space for the others that need to be examined.  
- The second type is the bam coverage plots and each file represents one quality level (file name start with bamcoverage);  
- The third type contains information that shows additional QC measures and here are the details:

  STAR_minimal_counts_soft_threshold: Total number of genes recovered was plotted as a function of STAR counts, and saturation function was used to fit the data. Once fitted, the algorithm first detected the minimal STAR counts required for yielding defined gene recovery percentage and compared it to the set level. The higher value of them was then used as the final minimal STAR counts threshold for separating samples that need to be resequenced (if less than the threshold and have no other issues). After that, normal distribution parameters were calculated for total gene numbers of all other samples that are above the minimal STAR counts threshold, and we considered any sample fell below 95% confidence interval need to be manually QCed. The minimal percentage of gene recovery can be specified in the json config file.

  scenario 1: a soft threshold was used
  <img src='img/Saturation_case1.png' width = 800>
  scenario 2: soft threshold less than set level and the later one was used
  <img src='img/Saturation_case2.png' width = 800>

    Spearman_correlation:Pairwise spearman correlation of all samples was calculated for detecting outlier. We assume the mean spearman value of each sample should be within the normal distribution of all sample means with 95% confidence (one tail). If not then that sample will be marked as 'outlier'.   

  <img src='img/SP_corr.png' width = 500>  
