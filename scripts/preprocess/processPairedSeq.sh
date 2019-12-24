#!/bin/sh

#BASEDIR=$(dirname "$0")
#echo "$BASEDIR"

#. $BASEDIR/getProperties.sh
#PROPERTY_FILE=$BASEDIR/genome.properties

echo "--- process sequence: Parse command line arguments"
while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
                -g|--genome_type)
                GENOME_TYPE="$2"
                shift # past argument
                ;;
	        -s|--pipeline_dir)
                PIPELINE_DIR="$2"
                shift # past argument
                ;;
                -i|--input_dir)
                INPUT_DIR="$2"
                shift # past argument
                ;;
                -o|--output_dir)
                OUTPUT_DIR="$2"
                shift # past argument
                ;;
	            -f|--fastqs)
		        FASTQS="$2"
                shift # past argument
                ;;
                *)
                # unknown option
                echo "Unknown option: $key, exiting."
                echo "Usage: ./processedPairedSeq.sh -g <genome_type> -s <pipeline_dir> -i <input_dir> -w <output_dir> -f <file array>"
                exit
                ;;
        esac
        shift # past argument or value
done

# get configuration values
PROPERTY_FILE=$PIPELINE_DIR/config/genome.properties
. $PIPELINE_DIR/scripts/getProperties.sh


HUMAN_GENOME=$(getProperty HUMAN_GENOME $PROPERTY_FILE)
HUMAN_STAR_INDEX=$(getProperty HUMAN_STAR_INDEX $PROPERTY_FILE)
HUMAN_GENOME_GTF=$(getProperty HUMAN_GENOME_GTF $PROPERTY_FILE)

MOUSE_GENOME=$(getProperty MOUSE_GENOME $PROPERTY_FILE)
MOUSE_STAR_INDEX=$(getProperty MOUSE_STAR_INDEX $PROPERTY_FILE)
MOUSE_GENOME_GTF=$(getProperty MOUSE_GENOME_GTF $PROPERTY_FILE)

trimmomatic_path=$(getProperty trimmomatic_path $PROPERTY_FILE)
trimmomatic_jar=$(getProperty trimmomatic_jar $PROPERTY_FILE)
pe_adpter=$(getProperty trimmomatic_pe_adpter $PROPERTY_FILE)
case "$GENOME_TYPE" in
 Human )
        GENOME_DIR=$HUMAN_STAR_INDEX
        GENOME_GTF=$HUMAN_GENOME_GTF
        echo $GENOME_DIR
        echo $GENOME_GTF;;
 *)
        GENOME_DIR=$MOUSE_STAR_INDEX
        GENOME_GTF=$MOUSE_GENOME_GTF
        echo $GENOME_DIR
        echo $GENOME_GTF;;
esac

echo "fastqs="$FASTQS

fastq_list=$(<$FASTQS)
IFS=',' read -ra FASTQ_ARRAY<<<"$fastq_list"


##split paired fastqs
cd $INPUT_DIR
echo "input="$INPUT_DIR
for file in ${FASTQ_ARRAY[@]}; do
    run_num=${file%.fastq}
    echo "runnum="$run_num
    fastq-dump --split-3 $run_num

    mv *_1.fastq $OUTPUT_DIR/paired_fastqs/
    mv *_2.fastq $OUTPUT_DIR/paired_fastqs/
done

#############Step 1: Perform QC and Trim Reads
cd $OUTPUT_DIR/paired_fastqs
for file in ${FASTQ_ARRAY[@]}; do
    fq=${file%.fastq}
    echo "fq=$fq"
    fq1_suffix="_1.fastq"
    fq2_suffix="_2.fastq"
    fq1=$fq$fq1_suffix
    fq2=$fq$fq2_suffix
    echo "fq1=$fq1"
    echo "fq2=$fq2"

    echo "Performing FastQC for original $fq"
    fastqc $fq1 -o $INPUT_DIR/fastQC_output
    fastqc $fq2 -o $INPUT_DIR/fastQC_output

    echo "Start to trim reads for $fq"
    java -jar $trimmomatic_path/$trimmomatic_jar PE -phred33 $fq1 $fq2 $fq"_1_p.fastq" $fq"_1_u.fastq" $fq"_2_p.fastq" $fq"_2_u.fastq" \
ILLUMINACLIP:$trimmomatic_path/$pe_adpter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

mv *_p.* $OUTPUT_DIR/trimmed_p
mv *_u.* $OUTPUT_DIR/trimmed_u

echo  "Started to QC Paired Reads"
cd $OUTPUT_DIR/trimmed_p
for file in ${FASTQ_ARRAY[@]}; do
    fq=${file%.fastq}
    echo $fq
    fq1_suffix="_1_p.fastq"
    fq2_suffix="_2_p.fastq"
    fq1=$fq$fq1_suffix
    fq2=$fq$fq2_suffix
    echo "Performing FastQC for trimmed $fq"
    fastqc $fq1 -o $OUTPUT_DIR/fastQC_output
    fastqc $fq2 -o $OUTPUT_DIR/fastQC_output
done
############Step 2: Align and assemble paired-end sequencing reads ##########
echo  "Started to Align Reads"
for file in ${FASTQ_ARRAY[@]}; do
    fq=${file%.fastq}
    echo $fq
    fq1_suffix="_1_p.fastq"
    fq2_suffix="_2_p.fastq"
    fq1=$fq$fq1_suffix
    fq2=$fq$fq2_suffix

    echo "Start to run mapping for $fq"
    STAR --runThreadN 16 --genomeDir $GENOME_DIR --sjdbGTFfile $GENOME_GTF --sjdbOverhang 100 --genomeSAsparseD 2 --readFilesIn $fq1 $fq2 --outFileNamePrefix ../star_output/$fq --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3

    sam_suffix="Aligned.out.sam"
    bam_suffix="Aligned.sorted.bam"
    in_name="$OUTPUT_DIR/star_output/$fq$sam_suffix"
    out_name=$fq$bam_suffix
    echo "Started to Sort Reads"
    samtools sort -n $in_name -o $out_name
done

mv *.bam $OUTPUT_DIR/star_output

#############Step 3: Feature Quantification##############
echo  "Started Feature Quantification"
cd $OUTPUT_DIR/star_output
for file in ${FASTQ_ARRAY[@]}; do
    echo $file
    fq=${file%.fastq}
    bam_suffix="Aligned.sorted.bam"
    bam_filename=$fq$bam_suffix
    ante="htseq_"
    suffix="_results.txt"
    out=$ante$fq$suffix
    echo $out
    htseq-count -f bam -s no $bam_filename $HUMAN_GENOME_GTF > $out
    mv $out $OUTPUT_DIR/htseqCount_output
done
