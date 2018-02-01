# TayWhale - a RNA seq pipeline
The TayWhale pipeline performs alignment using STAR, differetial expression analysis using Cufflinks, Fusion transcript detection using STAR-Fusion, and allele specific expression using GATK

The pipeline is still being developed!

# Command line
Activate the conda environment

        source activate TayWhale

Next, you may run the pipeline

    nextflow TayWhale.nf --r1 read1.fq --r2 --read2.fq --sample sampleID --output output_directory -c config

# Install
First install the conda  environment:

         ./create_conda_env.sh

NOTE: you need to install conda before running the script!

The conda script will create an environment named TayWhale. If you cannot use conda to install packages, you will have to install the dependencies manually. The dependencies include:

    STAR
    
    STAR-Fusion

    Samtools

    Picard tools

    CuffLinks

    RepeatMasker

    Blast

Then activate the conda environment, and run the install script:

        source activate TayWhale

        python install.py <reference.fasta> <GTF_file> <output_dir> <max_read_length> > Taywhale.conf

Note the script does not apply repeatmasking, if you want to mask repeats, you could either apply the repeatmasker to the reference genome before running the script. You need to suply the absolute path of your output directory.

The install script creates a folder named <output_dir> containg the indexed reference file and resources needed by the pipeline. The script also generates a config file (TayWhale.conf) this file may be edited to make the pipeline run on
slurm etc. visit the nextflow website for more info on how to setup the config:

    https://www.nextflow.io/docs/latest/config.html

Lastly, enter the config file and set the path to GATK. Example:

        nano TayWhale.conf

And edit:    
    
        GATK="/sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar"
