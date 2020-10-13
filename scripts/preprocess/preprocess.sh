##$ -S /bin/bash
##$ -cwd -V
#!/bin/sh
#default

echo "--- preprocess: Parse command line arguments"
while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
                 -g|--genome_type)
                GENOME_TYPE="$2"
                shift # past argument
                ;;
                -p|--is_paired)
                IS_PAIRED="$2"
                shift # past argument
                ;;
                -d|--is_dropseq)
                IS_DROPSEQ="$2"
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
	            -h|--help)
                echo "Usage: ./preprocess.sh -g <genome_type:Human or Mouse> -p <is paired end read: true or false>  -d <use Drop-Seq:true or false> -s <pipeline directory> \
                -i <fastq reads input directory> -o <preprocess output directory> -f <fastq list file>"
                exit
                ;;
                *)
                # unknown option
                echo "Unknown option: $key, exiting."
                echo "Usage: ./preprocess.sh -g <genome_type:Human or Mouse> -p <paired end reads: true or false>  -d <use Drop-Seq method:true or false>  -s <pipeline directory> \
                -i <fastq reads input directory> -o <preprocess output directory> -f <fastq list file>"
                exit
                ;;
        esac
        shift # past argument or value
done

SCRIPT_DIR=$PIPELINE_DIR/scripts
echo $SCRIPT_DIR

if [[ ! -d  $OUTPUT_DIR ]]; then
        echo "Creatig Output Directory: $OUTPUT_DIR"
        mkdir -p $OUTPUT_DIR
fi

cd $INPUT_DIR
#mkdir for store original fastq QC results
mkdir -p fastQC_output

cd $OUTPUT_DIR 
## create dirs if not exists
mkdir -p fastQC_output
mkdir -p star_output
mkdir -p htseqCount_output
mkdir -p trimmed_p
mkdir -p trimmed_u
#mkdir -p analysis_scripts
mkdir -p paired_fastqs

echo "is_drop="$IS_DROPSEQ
##############Analysis SRA Data ###############
cd $SCRIPT_DIR
if [ "$IS_DROPSEQ" = "true" ]; then
    ./processDropSeq.sh -g $GENOME_TYPE -s $PIPELINE_DIR -i $INPUT_DIR -o $OUTPUT_DIR -f $FASTQS
elif [ "$IS_PAIRED" = "true" ]; then
    ./processPairedSeq.sh -g $GENOME_TYPE -s $PIPELINE_DIR -i $INPUT_DIR -o $OUTPUT_DIR -f $FASTQS
else
    ./processSingleSeq.sh -g $GENOME_TYPE -s $PIPELINE_DIR -i $INPUT_DIR -o $OUTPUT_DIR -f $FASTQS
fi
