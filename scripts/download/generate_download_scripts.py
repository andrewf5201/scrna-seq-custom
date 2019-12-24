import os
import glob
import argparse
import sys

from string import Template

parser = argparse.ArgumentParser(description="generating raw data download scripts")

parser.add_argument('--inputdir', help='specify a full path of metadata file directory')
parser.add_argument('--prefix', help="specify metadata file name prefix")
parser.add_argument('--method', help="specify download method: wget or fastq-dump")
parser.add_argument('--destdir', help='specify a full path of download destination directory')

arguments = parser.parse_args()

dump_script_template = Template('$SCRIPT_DIR/download_with_fastq_dump.sh -m $metadata -d $destdir')
wget_script_template = Template('$SCRIPT_DIR/download_with_wget.sh -m $metadata -d $destdir')

if __name__ == "__main__":

    inputdir = arguments.inputdir
    prefix = arguments.prefix
    meta_files = inputdir+"/" +prefix +"_metadata_*"
    method = arguments.method
    destdir = arguments.destdir

    pipeline_dir = inputdir[0:inputdir.rfind('/') ]
    scriptdir = pipeline_dir+"/scripts"

    print (pipeline_dir)
    jobname='runDownloadJob.sh'
    job=open(jobname,'w')
    job.write('#!/bin/bash\n\n')

    for file in glob.glob(meta_files):
        start = file.rfind('_')
        end = file.rfind('.txt')
        suffix = file[start+1:end]
        output_filename = "download_"+prefix+"_"+suffix+".sh"
        print(output_filename)

        f=open(output_filename,"w")
        f.write('#!/bin/bash\n\n')

        if (method == 'wget'):
            f.write(wget_script_template.safe_substitute(SCRIPT_DIR=scriptdir, metadata=file, destdir=destdir))
        elif (method == 'fastq-dump'):
             f.write(dump_script_template.safe_substitute(SCRIPT_DIR=scriptdir, metadata=file, destdir=destdir))
        else:
            sys.exit('ERROR: Unknown download method')

        f.close()
        os.chmod(output_filename, 0o775)
        job.write('nohup bash '+ output_filename +' > ' + output_filename.replace('.sh', '.out') +' &\n')
    os.chmod(jobname, 0o775)
    job.close()
