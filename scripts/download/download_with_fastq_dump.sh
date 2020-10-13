#$ -S /bin/sh

while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
                 -m|--metadata-file)
                metadata="$2"
                shift # past argument
                ;;
                -d|--dest_dir)
                downloaddir="$2"
                shift # past argument
                ;;
                -p|--is_paired)
                IS_PAIRED="$2"
                shift # past argument
                ;;
                -h|--help)
                echo "Usage: ./download_with_fastq_dump.sh -m <metadata file full path> -d <download destination dir>  -p <is_paired>"
                exit
                ;;
                *)
                # unknown option
                echo "Unknown option: $key, exiting."
                exit
                ;;
        esac
        shift # past argument or value
done
echo $downloaddir

{
read
while IFS=$'\t' read f1 f2 f3 f4 f5
do
#  echo "1 : $f1"
#  echo "2 : $f2"
#  echo "3 : $f3"
#  echo "4 : $f4"
#  echo "5 : $f5"

fastqsdir=$downloaddir/$f2
mkdir -p $fastqsdir
cd $fastqsdir

echo "curr dir $fastqsdir"
RUN=$f4
echo "run=$RUN"
suffix=".fastq"
filename=$RUN$suffix
echo "fastq file="$filename

if [ -f "$filename" ]
then
  echo "skip downloading $name"
else
    if [ "$IS_PAIRED" = "true" ]; then
        fastq-dump --split-3 $RUN
    else
        fastq-dump $RUN
    fi
fi
cd $downloaddir
done } < $metadata
