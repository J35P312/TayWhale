import sys
import  os

programDirectory = os.path.dirname(os.path.abspath(__file__))
os.system("mkdir {}".format(sys.argv[3]))
os.system("cp {} {}/genome.fa".format(sys.argv[1],sys.argv[3]))
os.system("cp {} {}/transcripts.GTF".format(sys.argv[2],sys.argv[3]))

os.system("STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {}/genome.fa --sjdbGTFfile {} --sjdbOverhang {}".format(sys.argv[3],sys.argv[3],sys.argv[2],sys.argv[4]))

os.system("samtools faidx {}/genome.fa".format(sys.argv[3]))

os.system("picard CreateSequenceDictionary R=$3/genome.fa".format(sys.argv[3]))

os.system("wget https://github.com/FusionFilter/FusionFilter/archive/FusionFilterv0.4.0.zip")
os.system("unzip FusionFilterv0.4.0.zip")
os.system("./FusionFilter-FusionFilterv0.4.0/util/gtf_file_to_feature_seqs.pl  {}/transcripts.GTF {}/genome.fa cDNA > {}/cDNA_seqs.fa".format(sys.argv[3],sys.argv[3],sys.argv[3]))
os.system("makeblastdb -in {}/cDNA_seqs.fa -dbtype nucl".format(sys.argv[3]))

os.system("blastn -query {}/cDNA_seqs.fa -db {}/cDNA_seqs.fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-3 -lcase_masking -num_threads 16 -word_size 11  >  {}/blast_pairs.outfmt6".format(sys.argv[3],sys.argv[3],sys.argv[3]))
os.system("./FusionFilter-FusionFilterv0.4.0/prep_genome_lib.pl --genome_fa {}/genome.fa --gtf {}/transcripts.GTF --blast_pairs {}/blast_pairs.gene_syms.outfmt6.gz --cpu 16 ---output_dir {}/".format(sys.argv[3],sys.argv[3],sys.argv[3],sys.argv[3]))

template=""
for line in open("template_conf/FindSV.config"):
    template += line

template=template.replace("ref=\"\"","ref={}/genome.fa".format(sys.argv[3]))

template=template.replace("STAR_ref_dir=\"\"","STAR_ref_dir={}".format(sys.argv[3]))

template=template.replace("ctat_folder=\"\"","ctat_folder={}".format(sys.argv[3]))

template=template.replace("BootstrapAnn=","BootstrapAnn={}/BootstrapAnn/BootstrapAnn.py".format(programDirectory))

print template



