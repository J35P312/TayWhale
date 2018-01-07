mkdir $3
cp $1 $3/genome.fa
cp $2 $3/transcripts.GTF

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $3 --genomeFastaFiles $3/genome.fa --sjdbGTFfile $2 --sjdbOverhang $4

samtools faidx $3/genome.fa
picard CreateSequenceDictionary R=$3/genome.fa

wget https://github.com/FusionFilter/FusionFilter/archive/FusionFilterv0.4.0.zip
unzip FusionFilterv0.4.0.zip
./FusionFilter-FusionFilterv0.4.0/util/gtf_file_to_feature_seqs.pl  $3/transcripts.GTF $3/genome.fa cDNA > $3/cDNA_seqs.fa

#make the cDNA_seqs.fa file blastable
makeblastdb -in $3/cDNA_seqs.fa -dbtype nucl
# perform the blastn search
blastn -query $3/cDNA_seqs.fa -db $3/cDNA_seqs.fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-3 -lcase_masking -num_threads 16 -word_size 11  >  $3/blast_pairs.outfmt6

#generate ctat folder
prep_genome_lib.pl --genome_fa $3/genome.fa --gtf $3/transcripts.GTF --blast_pairs $3/blast_pairs.gene_syms.outfmt6.gz --cpu 16 ---output_dir $3/

#copy the config template
cp template_conf/FindSV.config TayWhale.conf

#replace the reference variables
ref="ref=\"\""
path="ref=\"$3\/genome.fa\""
sed -i -e "s/$ref/$path/g" TayWhale.conf

ref="ref=\"\""
path="ref=\"$3\/genome.fa\""
sed -i -e "s/$ref/$path/g" TayWhale.conf

ref="STAR_ref_dir=\"\""
path="STAR_ref_dir=\"$3\""
sed -i -e "s/$ref/$path/g" TayWhale.conf

ref="ctat_folder=\"\""
path="ctat_folder=\"$3\""
sed -i -e "s/$ref/$path/g" TayWhale.conf
