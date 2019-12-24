import configparser
import os
import getopt
import sys
from string import Template

def usage():
    print ("Usage: " + sys.argv[0] +" --config   config_file_name\n")

    print ("Example:")
    print ("\t" +  sys.argv[0] +" --config ../../config/config.properties")
    sys.exit(1)

def read_properties_file(prop_file):
    with open(prop_file, 'r') as f:
        config_string = '[config_section]\n' + f.read()
        config = configparser.ConfigParser()
        config.read_string(config_string)
        return dict(config.items('config_section'))

if __name__ == "__main__":
    config_file = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:h", ["config="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    for o, a in opts:
     if o in ("-c", "--config"):
         config_file = a
     elif o in ("-h", "--help"):
         usage();
         sys.exit()
     else:
         assert False, "Unknown option"
         sys.exit(-1)

    props = read_properties_file(config_file)

    genome_type = props['genome_type']
    is_paired = props['is_paired']
    fastqs_input_dir = props['fastqs_input_dir']
    pipeline_dir = props['pipeline_dir']
    experiment = props['experiment']
    preprocess_dir = pipeline_dir +'/preprocess/' +experiment
    script_dir = pipeline_dir +'/scripts'
    config_dir = pipeline_dir +'/config'

    batch_size_str = props.get("batch_size")

    if batch_size_str is None:
        batch_size = 100 #default number
    else:
        batch_size = int(props['batch_size'])


    print(preprocess_dir)
    #dir_path = os.path.dirname(os.path.realpath(__file__))

    individual_dirs = os.listdir(fastqs_input_dir)
    individual_dirs.sort()

    script_template = Template('$SCRIPT_DIR/preprocess/preprocess.sh -g $genome_type -p $is_paired -s $PIPELINE_DIR -i '
                               '$FASTQS_INPUT_DIR/$individual -o $PREPROCESS_DIR/$individual '
                               '-f $FASTQS_INPUT_DIR/$individual/$fastqs')
    multiqc_template = Template('multiqc $FASTQS_INPUT_DIR/$individual/fastQC_output $FASTQS_INPUT_DIR/$individual/star_output'
                                 ' $FASTQS_INPUT_DIR/$individual/htseqCount_output --filename $individual_multiqc_report.html')
    count=0
    jobname='runPreprocessJob.sh'
    job=open(jobname,'w')
    job.write('#!/bin/bash\n\n')

    multiqc_jobname="runMultiQCJob.sh"
    multiqc_job=open(multiqc_jobname,'w')
    multiqc_job.write('#!/bin/bash\n\n')

    for individual in individual_dirs:
        if os.path.isdir(os.path.join(fastqs_input_dir,individual)):
            print("individual="+individual)
            multiqc_filename='runMultiqc_'+individual +'.sh'
            multiqc_file=open(multiqc_filename, "w")
            multiqc_file.write('#!/bin/bash\n\n')
            multiqc_file.write(multiqc_template.safe_substitute(FASTQS_INPUT_DIR=fastqs_input_dir, individual=individual))
            multiqc_file.close()
            multiqc_job.write('SBATCH ' + multiqc_filename + '\n')

            fastqarr=[]
            #init fastq array
            for file in os.listdir(os.path.join(fastqs_input_dir,individual)):
                if file.endswith(".fastq"):
                    fastqarr.append(file)
                    fastqarr.sort()

            #split fastq array based on the batch size
            i=0
            print("fastqarr="+str(fastqarr))
            while i<len(fastqarr):
                partition=fastqarr[i:i+batch_size]
                file_seq = str(int(i/batch_size+1))
                fastqs_filename=individual+"_"+file_seq+".txt"
                fastqs_file=open(os.path.join(fastqs_input_dir,individual,fastqs_filename),'w')
                print(fastqs_filename)
                print(",".join(partition))
                fastqs_file.write(",".join(partition))
                fastqs_file.close()

                #generate preprocess sh file for each partition
                filename='preprocess_'+individual+'_'+file_seq+'.sh'
                f=open(filename,"w")
                f.write('#!/bin/bash\n\n')
                f.write(script_template.safe_substitute(SCRIPT_DIR=script_dir, genome_type=genome_type,
                                                        is_paired=is_paired, PIPELINE_DIR=pipeline_dir,
                                                        FASTQS_INPUT_DIR=fastqs_input_dir, individual=individual,
                                                        PREPROCESS_DIR=preprocess_dir, fastqs=fastqs_filename))
                f.close()
                os.chmod(filename, 0o775)
                i += batch_size
                job.write('SBATCH '+ filename + "\n")
    os.chmod(jobname, 0o775)
    os.chmod(multiqc_jobname, 0o775)
    job.close()
    multiqc_job.close()

