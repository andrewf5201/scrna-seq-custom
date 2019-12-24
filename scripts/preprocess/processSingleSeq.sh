#$ -S /bin/sh
# -cwd -V

echo "--- process  sequence: Parse command line arguments"
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
                echo "Usage: ./processedSingleSeq.sh -g <genome_type> -s <script_dir> -i <input_dir> -o <output_dir> -f <fastq list file>"
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
se_adpter=$(getProperty trimmomatic_se_adpter $PROPERTY_FILE)


#read comma delimited fastq list
fastq_list=$(<$FASTQS)
IFS=',' read -ra FASTQ_ARRAY<<<"$fastq_list"
echo "fastq_list="$fastq_list

shopt -s nocasematch
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

#############Step 1: Perform Trim Reads and QC
echo  "------Started to QC and  Trim Reads---------"
cd $INPUT_DIR
shopt -s nullglob
for file in ${FASTQ_ARRAY[@]}; do
  echo "Performing FastQC for original $file"
  fastqc $file -o $INPUT_DIR/fastQC_output

  fq=${file%.fastq}
  echo "runnum="$fq
  outfile=$fq"_u.fastq"
  echo "Start to trim reads for $fq"
  java -jar $trimmomatic_path/$trimmomatic_jar SE -phred33 $file $OUTPUT_DIR/trimmed_u/$outfile ILLUMINACLIP:$trimmomatic_path/$se_adpter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

trimmed_array=()
for i in ${FASTQ_ARRAY[@]};do
    base=$(echo $i| cut -f1 -d '.')
    trim_suffix="_u.fastq"
    trimmed_array+=($base$trim_suffix)
done

echo  "------Started to QC Reads after Trimming--------"
cd $OUTPUT_DIR/trimmed_u
for fq in ${trimmed_array[@]}; do
    echo "Performing FastQC for $fq"
    fastqc $fq -o $OUTPUT_DIR/fastQC_output
done
############Step 2: Align and assemble single-end sequencing reads ##########
echo  "--------Started to Align Reads----------"
for fq in ${trimmed_array[@]}; do
    fqbase=$(echo $fq| cut -f1 -d '_')
    #echo $fqbase
    echo "Start to run mapping for $fq"
    STAR --runThreadN 1 --genomeDir $GENOME_DIR --sjdbGTFfile $GENOME_GTF --sjdbOverhang 100 --genomeSAsparseD 2 --readFilesIn $fq --outFileNamePrefix ../star_output/$fqbase --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3

    sam_suffix="Aligned.out.sam"
    bam_suffix="Aligned.sorted.bam"
    in_name="$OUTPUT_DIR/star_output/$fqbase$sam_suffix"
    out_name=$fqbase$bam_suffix
    echo "Started to Sort Reads"
    samtools sort -n $in_name -o $out_name
done

mv *.bam $OUTPUT_DIR/star_output

#############Step 3: Feature Quantification##############
STAR_ARRAY=()
for i in ${FASTQ_ARRAY[@]}; do
    base=$(echo $i| cut -f1 -d '.')
    bam_suffix="Aligned.sorted.bam"
    STAR_ARRAY+=($base$bam_suffix)
done

echo  "------Started Feature Quantification----------"
cd $OUTPUT_DIR/star_output
for file in ${STAR_ARRAY[@]}; do
    echo $file
    prefix=$(echo $file | cut -f1 -d "A" )
    echo $prefix
    ante="htseq_"
    suffix="_results.txt"
    out=$ante$prefix$suffix
    echo $out
    htseq-count -f bam -s no $file $HUMAN_GENOME_GTF > $out
    mv $out $OUTPUT_DIR/htseqCount_output
done
