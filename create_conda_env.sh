conda config --add channels r
conda config --add channels bioconda

conda create --force -n TayWhale python=2.7 scipy numpy
source activate TayWhale

conda install -c bioconda STAR --yes
conda install -c bioconda picard --yes
conda install -c bioconda samtools --yes
conda install -c bioconda STAR-Fusion --yes
conda install -c bioconda cufflinks --yes
conda install -c bioconda trinity --yes
conda install -c bioconda salmon --yes

source deactivate
