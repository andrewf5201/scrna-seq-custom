#$ -S bin/sh
# -cwd -V

while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
                -p|--pipeline_dir)
                pipeline_dir="$2"
                shift # past argument
                ;;
                -c|--config_file)
                config_file="$2"
                shift # past argument
                ;;
                -h|--help)
                echo "Usage: ./runAnalysis.sh -p <pipeline directory> -c <config file>"
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

#reading properties
. $pipeline_dir/scripts/getProperties.sh

CONF_PROPERTY_FILE=$pipeline_dir/config/$config_file
experiment=$(getProperty "experiment"  $CONF_PROPERTY_FILE)
metadata_file=$(getProperty "metadata_file" $CONF_PROPERTY_FILE)
count_matrix=$(getProperty "matrix_file" $CONF_PROPERTY_FILE)
build_count_matrix=$(getProperty "build_count_matrix" $CONF_PROPERTY_FILE)

echo "matrix="$count_matrix
echo $pipeline_dir
cd $pipeline_dir


#absolute path
SCRIPT_DIR=$pipeline_dir/scripts/analysis
PREPROCESS_DIR=$pipeline_dir/preprocess/$experiment
ANALYSIS_DIR=$pipeline_dir/analysis/$experiment

#echo $PREPROCESS_DIR
if $build_count_matrix; then
  echo "build count matrix file"
  if [ ! -f "$count_matrix" ]; then
      for d in  $PREPROCESS_DIR/*
      do
         if [ -d $d ]
         then
           individual=${d##*/}
           echo $individual
           echo "subdir="$d
           echo ${individual}_counts.txt
           #create count matrix for each individual
           Rscript $SCRIPT_DIR/createCountFile.R $PREPROCESS_DIR/$individual/htseqCount_output $ANALYSIS_DIR/data ${individual}_counts.txt $metadata_file
         fi
      done
      #assemble all individual count file
      Rscript $SCRIPT_DIR/buildCountMatrix.R $ANALYSIS_DIR/data $count_matrix
  fi
fi



#setup analysis directory
mkdir -p analysis
mkdir -p analysis/$experiment
mkdir -p analysis/$experiment/data

resultdir=analysis/$experiment/results
mkdir -p $resultdir
mkdir -p $resultdir/scater
mkdir -p $resultdir/edgeR
mkdir -p $resultdir/DESeq

norm_methods=$(getProperty "norm_methods"  $CONF_PROPERTY_FILE)
echo $norm_methods

IFS=',' read -ra array <<< "$norm_methods"

for method in "${array[@]}"
do
    echo "$method"
    packageName=${method%:*}
    algorithm=${method#*:}
    echo "package="$packageName
    echo "method="$algorithm

    #to lower case
    packageName="$(tr [A-Z] [a-z] <<< "$packageName")"
    algorithm="$(tr [A-Z] [a-z] <<< "$method")"


    if [ "$packageName" = "deseq" ];
    then
        echo DESeq
        Rscript $SCRIPT_DIR/Deseq.R  $algorithm $CONF_PROPERTY_FILE
    elif [ "$packageName" = "scater" ]
    then
        echo scater
        Rscript $SCRIPT_DIR/scater_2.R  $algorithm $CONF_PROPERTY_FILE
    elif [ "$packageName" = "edger" ]
    then
        echo edgeR
        Rscript $SCRIPT_DIR/edgeR.R  $algorithm $CONF_PROPERTY_FILE
    else
        echo "unknown package name:$packageName, valid package names are DeSEQ, scater, and edgeR."
    fi
done
