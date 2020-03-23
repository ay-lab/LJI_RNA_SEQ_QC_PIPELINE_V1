# **Reference genome**

---

The general rule of making reference index should follow [STAR's tutorial](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), and any change/addition/deletion of the reference files must be properly documented.

Here are the source of implemented reference files where we downloaded reference genome from:
1. Human genome version GRCh37: **Full** genome reference downloaded from [GENCODE Release 19 (GRCh37.p13)](https://www.gencodegenes.org/human/release_19.html), Bed file downloaded from [Rseqc reference hg19](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_GencodeCompV19.bed.gz/download).

2. Human genome version GRCh38: **Primary assemble** genome reference downloaded from [GENCODE Release 32 (GRCh38.p13)](https://www.gencodegenes.org/human/release_32.html), Bed file downloaded from [Rseqc reference hg38](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_Gencode_V28.bed.gz/download).

3. Mouse genome version GRCm38:  **Primary assemble** genome reference downloaded from [GENCODE Release 23](https://www.gencodegenes.org/mouse/release_M23.html), Bed file downloaded from [Rseqc reference mm10](https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_Gencode_VM18.bed.gz/download).


Importantly, STAR recommended to exclude patches and alternative haplotypes from the genome model (see tutorial section 2.2.1), therefore you should only use the primary assemble for mapping. In case only full genome model is available, you can use the following python script or your own code to trim the genome file.

- touch remove_patch.py and put in the following code

```
# /usr/bin/env python3

# remove PATCH and Alternative haplotypes from fasta file
# usage remove_patch.py in.fasta > output.fasta


from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
ffile = SeqIO.parse(in_fasta, "fasta")
header_pattern = ['PATCH','HSCHR']
for seq_record in ffile:
    if not any([i in seq_record.description for i in header_pattern]):
        print(seq_record.format('fasta'))
```



- For making reference index; please make proper changes according to your settings

```
#!/bin/bash
#PBS -N STAR_gen_37
#PBS -o /mnt/BioScratch/ndu/gen_reference/out_STAR_gen_37
#PBS -e /mnt/BioScratch/ndu/gen_reference/err_STAR_gen_37
#PBS -q default
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb
#PBS -l walltime=20:00:00
cd  /mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference
mkdir -p /mnt/BioScratch/ndu/gen_reference

# trim off batches and alternative haplotypes if needed
./remove_patch.py GRCH37.P13/GRCh37.p13.genome.fa > GRCH37.P13/GRCh37.p13.genome.primary_assembly.fa

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./GRCH37.P13 --genomeFastaFiles ./GRCH37.P13/GRCh37.p13.genome.primary_assembly.fa --sjdbGTFfile ./GRCH37.P13/gencode.v19.annotation.gtf --sjdbOverhang 100
```

Separately, an annotation file should be made for counting reads by gene type (gene_type_4) and TPM calculation in the pipeline. A example human GRCh37 annotation file can be downloaded [here](./files/GRCh37_annotation.csv). To make the annotation table, you will need to execute the following steps:

1. Getting gene name and type for each ensembl ID. You can export annotations table for [GRCh37](https://grch37.ensembl.org/biomart/martview/923bfa0c7a1727c3fe634eb8c422df78) and [GRCh38](http://uswest.ensembl.org/biomart/martview/e859cf18a85550949a12ba09c8ab117c) from ensembl biomart. Mouse gene annotations are also available from these links.
2. Merge gene types so 4 categories. The dictionary for merging of gene types can be downloaded [here](files/dict_gene_type.csv) .
3. Getting exon length for each gene; this is a common output from most counting tools such as featureCounts and HTSeq, or you can download from a confident source.
