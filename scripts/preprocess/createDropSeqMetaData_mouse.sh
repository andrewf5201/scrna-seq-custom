#!/bin/bash 
#$ -cwd -V
#$ -S /bin/sh


DROP_SEQ_TOOLS_PATH=/home/andrew/tools/Drop-seq_tools-2.3.0
SAM_TOOLS_PATH=/opt/wangcluster/samtools/1.9/bin
BGZIP_PATH=/home/andrew/tools/htslib-1.5/bin

GENOME_HOME=/home/andrew/genome/GRCm38
reference_name=GRCm38.meta
reference_fasta=${GENOME_HOME}/GRCm38.primary_assembly.genome.fa
species=GRCm38
gtf=${GENOME_HOME}/gencode.vM15.primary_assembly.annotation.gtf
dropseq_root=$DROP_SEQ_TOOLS_PATH
output_dir=${GENOME_HOME}/meta
bgzip_executable=${BGZIP_PATH}/bgzip
samtools_executable=${SAM_TOOLS_PATH}/samtools


mkdir -p ${output_dir}

${DROP_SEQ_TOOLS_PATH}/create_Drop-seq_reference_metadata.sh -n ${reference_name} -r ${reference_fasta} -s ${species} -g ${gtf} -d ${dropseq_root} -o ${output_dir}  -b ${bgzip_executable} -i ${samtools_executable}
~                        
