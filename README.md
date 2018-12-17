# TayWhale - a RNA seq pipeline
The TayWhale pipeline performs alignment using STAR, transript quantification using Salmon, Fusion transcript detection using STAR-Fusion, transcript assembly using Sringie, and allele specific expression using GATK.

The pipeline is still being developed!

# Command line
The output of each tool will be stored in the <ouput_directory> (note, the fastqs must be gziped)

    nextflow TayWhale.nf --r1 read1.fq.gz --r2 --read2.fq.gz --sample sampleID --output output_directory -c config

# Install
Install the dependencies:

	singularity
	vep
	python 2.7, numpy, scikit
	GATK

Create indexes for:

	salmon
	star fusion
	star
	
remember to generate genome dictionary and fasta index files!

Download or create the singularity image:

	singularity pull --name customname.img shub://J35P312/TayWhale
	sudo singularity build <imagename.simg> SINGULARITY

run the install script to generate a config file:

	python install.py

Lastly, enter the config file and set the path variables (GATK; references, vep, and singularity image):

        nano TayWhale.conf


![TayWhale](TayWhale.jpg)
