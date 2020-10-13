##$ -S /bin/bash
##$ -cwd -V
#!/bin/sh
#default

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
                echo "Usage: ./processedDropSeq.sh -g <genome_type> -s <pipeline_dir> -i <input_dir> -w <output_dir> -f <file array>"
                exit
                ;;
        esac
        shift # past argument or value
done

# get configuration values
PROPERTY_FILE=$PIPELINE_DIR/config/genome.properties
. $PIPELINE_DIR/scripts/getProperties.sh


HUMAN_GENOME_META_DIR=$(getProperty HUMAN_GENOME_META_DIR $PROPERTY_FILE)
HUMAN_REFERENCE_FASTA=$(getProperty HUMAN_REFERENCE_FASTA $PROPERTY_FILE)

MOUSE_GENOME_META_DIR=$(getProperty MOUSE_GENOME_META_DIR $PROPERTY_FILE)
MOUSE_REFERENCE_FASTA=$(getProperty MOUSE_REFERENCE_FASTA $PROPERTY_FILE)

DROP_SEQ_PATH=$(getProperty DROP_SEQ_PATH $PROPERTY_FILE)
NUM_CORE_BARCODE=$(getProperty NUM_CORE_BARCODE $PROPERTY_FILE)
MIN_NUM_GENES_PER_CEL=$(getProperty MIN_NUM_GENES_PER_CEL $PROPERTY_FILE)

case "$GENOME_TYPE" in
 Human )
        GENOME_DIR=$HUMAN_GENOME_META_DIR
        GENOME_REF_FASTA=$HUMAN_REFERENCE_FASTA
        echo $GENOME_DIR
        echo $GENOME_REF_FASTA;;
 *)
        GENOME_DIR=$MOUSE_GENOME_META_DIR
        GENOME_REF_FASTA=$MOUSE_REFERENCE_FASTA
        echo $GENOME_DIR
        echo $GENOME_REF_FASTA;;
esac

echo "fastqs="$FASTQS

fastq_list=$(<$FASTQS)
IFS=',' read -ra FASTQ_ARRAY<<<"$fastq_list"

cd $INPUT_DIR
echo "input="$INPUT_DIR
for file in ${FASTQ_ARRAY[@]}; do
    echo $file
    run_num=${file%.fastq}
    echo "runnum="$run_num

    mv ${run_num}_1.fastq $OUTPUT_DIR/paired_fastqs/
    mv ${run_num}_2.fastq $OUTPUT_DIR/paired_fastqs/
done

# create output directories
mkdir -p $OUTPUT_DIR/bam
mkdir -p $OUTPUT_DIR/dge
mkdir -p $OUTPUT_DIR/read_counts

cd $OUTPUT_DIR/paired_fastqs
for file in ${FASTQ_ARRAY[@]}; do
    fq=${file%.fastq}
    echo $fq
    fq1_suffix="_1.fastq"
    fq2_suffix="_2.fastq"
    fq1=$fq$fq1_suffix
    fq2=$fq$fq2_suffix

    echo "Convert fastq to bam file"

    java -jar $DROP_SEQ_PATH/3rdParty/picard/picard.jar FastqToSam \
      F1=$fq1 \
      F2=$fq2 \
      O=$OUTPUT_DIR/bam/$fq.bam \
      SM=$fq

    mkdir -p $OUTPUT_DIR/align$fq
    mkdir -p $OUTPUT_DIR/tmp$fq

    echo "-----step1: Run Drop-seq_alignment.sh"
    echo $DROP_SEQ_PATH/Drop-seq_alignment.sh -g $GENOME_DIR/STAR -r $GENOME_REF_FASTA -d $DROP_SEQ_PATH -o $OUTPUT_DIR/align$fq -t $OUTPUT_DIR/tmp$fq  $OUTPUT_DIR/bam/$fq.bam
    $DROP_SEQ_PATH/Drop-seq_alignment.sh -g $GENOME_DIR/STAR -r $GENOME_REF_FASTA -d $DROP_SEQ_PATH -o $OUTPUT_DIR/align$fq -t $OUTPUT_DIR/tmp$fq  $OUTPUT_DIR/bam/$fq.bam

    echo "-----step2: Run BamTagHistogram"
    echo $DROP_SEQ_PATH/BamTagHistogram I=$OUTPUT_DIR/align$fq/final.bam O=$OUTPUT_DIR/read_counts/${fq}_readcounts.txt.gz TAG=XC FILTER_PCR_DUPLICATES=true
    $DROP_SEQ_PATH/BamTagHistogram I=$OUTPUT_DIR/align$fq/final.bam O=$OUTPUT_DIR/read_counts/${fq}_readcounts.txt.gz TAG=XC FILTER_PCR_DUPLICATES=true

    echo "----step3 Run DigitalExpression"
    echo $DROP_SEQ_PATH/DigitalExpression I=$OUTPUT_DIR/align$fq/final.bam O=$OUTPUT_DIR/dge/${fq}_dge.txt.gz SUMMARY=$OUTPUT_DIR/dge/${fq}_dge.summary.txt NUM_CORE_BARCODES=$NUM_CORE_BARCODE MIN_NUM_GENES_PER_CELL=$MIN_NUM_GENES_PER_CEL
    $DROP_SEQ_PATH/DigitalExpression I=$OUTPUT_DIR/align$fq/final.bam O=$OUTPUT_DIR/dge/${fq}_dge.txt.gz SUMMARY=$OUTPUT_DIR/dge/${fq}_dge.summary.txt NUM_CORE_BARCODES=$NUM_CORE_BARCODE MIN_NUM_GENES_PER_CELL=$MIN_NUM_GENES_PER_CEL
done
