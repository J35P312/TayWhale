# TayWhale - a RNA seq pipeline
The TayWhale pipeline performs alignment using STAR, differetial expression analysis using Cufflinks, Fusion transcript detection using STAR-Fusion, and allele specific expression using GATK

The pipeline is still bein developed!

# Command line
Activate the conda environment

        source ac

    nextflow TayWhale.nf --r1 read1.fq --r2 --read2.fq --sample sampleID --output output_directory --ref STAR_reference_folder -c config

# Install
Download STAR-Fusion
        
        https://github.com/STAR-Fusion/STAR-Fusion

run the generate_STAR_resource.sh script to create your reference  package (or download a premade via the star fusion Github). Type the  following to run the script:

       mkdir ref_folder

       cd ref_folder

       /path/generate_STAR_resource.sh genome.GTF genome.FA

The script will take about one day, and will consume 30GB harddrive space. The script requires repeatmasker, samtools, and blast to be installed. You must also set the STAR_FUSION_DIR parameter on line 16. This is the path of the STAR-fusion folder. Once the script is complete, you may use ref_folder as a STAR reference package. Note: the script performs the steps described here:
        
        https://github.com/FusionFilter/FusionFilter/wiki/Building-a-Custom-FusionFilter-Dataset

Visit this website if you run in to troubles creating your reference package.

    Setup The conda environment:

         ./create_conda_env.sh

    NOTE: you need to install conda before running the script!

The conda script will create an environment named TayWhale. If you cannot use conda to install packages, you will have to install the dependencies manually. The dependencies include:

    STAR
    
    STAR-Fusion

    Samtools

    Picard tools

    GATK

    CuffLinks
