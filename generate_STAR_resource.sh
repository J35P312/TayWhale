#!/bin/bash -l

#SBATCH -A sens2017106
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 74:00:00
#SBATCH -J STAR_index

#argument 1: the genome GTF
#argument 2: the reference genome

#howto: make a directory where you want to put your genome resources. cd into that directory, and run this script.
#./generate_STAR_resource.sh GTF genome.fa
#remember to change the STAR_FUSION_DIR_parameter!

STAR_FUSION_DIR=../STAR-Fusion-v1.2.0/

#get the exons of protein coding genes
grep -E "#|    exon    " $1 | grep -E "#|protein_coding" > transcripts.GTF
ln -s $2 genome.fa

#filter out the  sequence of transcripts
$STAR_FUSION_DIR/FusionFilter/util/gtf_file_to_feature_seqs.pl  transcripts.GTF genome.fa cDNA > cDNA_seqs.fa

#repeatmasking
RepeatMasker -pa 6 -s -species human -xsmall cDNA_seqs.fa

#make the cDNA_seqs.fa file blastable
makeblastdb -in cDNA_seqs.fa.masked -dbtype nucl
# perform the blastn search
blastn -query cDNA_seqs.fa.masked -db cDNA_seqs.fa.masked -max_target_seqs 10000 -outfmt 6 -evalue 1e-3 -lcase_masking -num_threads 16 -word_size 11  >  blast_pairs.outfmt6

#generate star index
$STAR_FUSION_DIR/FusionFilter/prep_genome_lib.pl --genome_fa genome.fa --gtf transcripts.GTF --blast_pairs blast_pairs.gene_syms.outfmt6.gz --cpu 16
