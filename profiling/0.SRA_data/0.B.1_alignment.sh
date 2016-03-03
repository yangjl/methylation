
#bowtie2/2.2.5

#http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
#Bismark

#(I) Running bismark_genome_preparation

module load bismark/0.14.3
module load bowtie2/2.2.5
bismark_genome_preparation --bowtie2 /home/jolyang/dbcenter/AGP/AGPv2/

#(II) Running bismark 

bismark --bowtie2 -n 1 -l 50 /data/genomes/homo_sapiens/GRCh37/ test_dataset.fastq


bismark --bowtie1 -n 1 /home/jolyang/dbcenter/AGP/AGPv2 -1 test_1.fastq -2 test_2.fastq

#(III) Running the Bismark bismark_methylation_extractor
bismark_methylation_extractor -s --comprehensive test_dataset.fastq_bismark.sam