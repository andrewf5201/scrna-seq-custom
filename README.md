# scrna-seq-custom
# Introduction
# Step 1: Define your Configuration File 
There are 2 configuration files the user needs to complete before the program can be run. *config.properties* defines certain aspects of each run, such as what analysis packages are used, or whether the reads used are single-ended or pair-ended. *genome.properties* defines the paths of the necessary genome and annotation files, both of which are needed for preprocessing.

The user needs to define the pipeline configuration properties in a “name=value” pair format.  

Sample Code for *config.properties*

```
#----- Dataset properties 
# specifies dataset’s genome type: human or mouse 
genome_type=human  
# indicates if the raw reads is paired-end or not 
is_paired=true  
# run accession prefix (SRR, ERR or other) 
read_prefix=SRR 
# user defined identifier for the experiment
experiment=BC 
 
#---- preprocessing  properties 
# specifies pipeline’s absolute path 
pipeline_dir=/home/pipeline 
# Raw read (fastq) files location 
fastqs_input_dir=/home/fastq_download/BC 
# Metadata file location 
metadata_file=/home/pipeline/metadata/GSE75688_metadata.txt 
# Number of fastq files to be downloaded and preprocessed in each batch 
batch_size=100 
 
#------Analysis Properties 
# indicates if need to build count matrix from preprocessing program generated count files build_count_matrix=true 
# defines count matrix file location. Used as the input for downstream analysis 
matrix_file=/home/pipeline/analysis/BC/data/BC_counts_matrix.txt 
# absolute path to a file contains human gene length data 
human_gene_length=/home/pipeline/config/hg38_gene_length.txt 
# a common delimited string contains one or many analysis methods in [package_name]:[algorightm] 
format norm_methods=DESeq:default 
# provides the column name of a factor in the metadata to perform analysis on 
specified_factor=cell_type 
```
Sample Code for *genome.properties*
```
# human reference genome file location 
HUMAN_GENOME=/home/genome/GRCh38 
# human genome STAR index file location 
HUMAN_STAR_INDEX=/home/genome/GRCh38/STAR_index/ 
# human Gencode annotation file location 
HUMAN_GENOME_GTF=/home/genome/GRCh38/gencode.v26.annotation.gtf 
 
#---Mouse genome properties 
MOUSE_GENOME=/home/genome/GRCm38 
MOUSE_STAR_INDEX=/home/genome/GRCm38/STAR_index/ 
MOUSE_GENOME_GTF=/home/genome/GRCm38/gencode.vM15.primary_assembly.annota tion.gtf 
 
#---Trimmomatic trimming tool properties  
trimmomatic_path=/home/Trimmomatic-0.36 trimmomatic_jar=trimmomatic-0.36.jar 
trimmomatic_pe_adpter=adapters/TruSeq3-PE.fa trimmomatic_se_adpter=adapters/TruSeq3-SE.fa 
```
# Step 2 (Optional): Downloading your Raw Read Files
If you do not already have the raw read files somewhere on your computer, and the files can be downloaded from an online database (e.g. Gene Expression Omnibus), You have the option of using the program itself to download the files for you.

