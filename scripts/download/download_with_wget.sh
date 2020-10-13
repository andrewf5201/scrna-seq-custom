#$ -S /bin/sh

echo "--- download: Parse command line arguments"
while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
                 -m|--metadata-path)
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
                echo "Usage: ./download_with_wget.sh -m <metadata file full path> -d <download destination dir> -p <is_paired>"
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
URL=$f5
URL=${URL%$'\r'}
name=${URL##*/}
echo "url=$URL"
echo "name=$name"
fastq_filename=${name%.gz}
if [ -f "$fastq_filename" ]
then
  echo "skip downloading $name"
else
  wget -c $URL
  if [[ "$name" == *.gz ]];
  then
      gunzip $name
      if [ -e $name ]
      then
          rm $name
      fi
   fi
fi

if [ "$IS_PAIRED" = "true" ]; then
     run_num=${fastq_filename%.fastq}
    ./deinterleave_fastq.sh deinterleave_paste.sh < $fastq_filename $run_num_1.fastq $run_num_2.fastq
fi
cd $downloaddir
done } < $metadata
