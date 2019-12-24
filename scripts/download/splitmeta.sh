#!/bin/sh
echo "---split metadata file: Parse command line arguments"
while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
            -d|--metafile_dir)
            metafile_dir="$2"
            shift # past argument
            ;;
            -f|--metadata_filename)
            metafile="$2"
            shift # past argument
            ;;
            -s|--batch_size)
            batchsize="$2"
            shift # past argument
            ;;
            -p|--output_prefix)
            prefix="$2"
            shift # past argument
            ;;
            -h|--help)
            echo "Usage: ./splitmeta.sh -d <metafile input directory> -f <metafile name> -s <batch size> -p<output filename prefix>"
            exit
            ;;
            *)
            # unknown option
            echo "Unknown option: $key, exiting."
            echo "Usage: ./splitmeta.sh -d <metafile input directory> -f <metafile name> -s <batch size> -p<output filename prefix>" 
            exit
            ;;
        esac
        shift # past argument or value
done
echo "input="$metafile_dir
echo "file="$metafile
echo "size="$batchsize

output_prefix=${prefix}_metadata_
echo $output_prefix
cd $metafile_dir
#tail -n +2 $metafile | split -d -l $batchsize - $output_prefix
tail -n +2 $metafile | split  -l $batchsize - $output_prefix

for file in ${output_prefix}*
do
    head -n 1 $metafile > tmp_file
    cat "$file" >> tmp_file
    mv -f tmp_file "$file"
done
