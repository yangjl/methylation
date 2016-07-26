### Jinliang Yang
### July 25th, 2016


library(farmeR)

idtab <- read.csv("data/teo20_ids.csv")
idtab$idchar <- gsub("-", "", idtab$idchar)
idtab$pid <- paste0("JR", idtab$plate)

### test: -i in the form of UCPAC ambiguity codes
"bcftools consensus -i -s NA001 -f in.fa in.vcf.gz > out.fa"

#JRIAL2A	JRIAL2B	JRIAL2C	JRIAL2D	JRIAL2E	JRIAL2F	JRIAL2G	JRIAL2H	JRIAL2I	JRIAL2J	
#JRIAL2K	JRIAL2L	JRIAL2M	JRIAL2N	JRIAL2O	JRIAL2P	JRIAL2Q	JRIAL2R	JRIAL2S	JRIAL2T

for(i in 1:nrow(idtab)){
    shid <- paste0("slurm-script/run_pg", i, ".sh")
    
    cmd1 <- paste0("mkdir largedata/wgbs_pgenome/", idtab$idchar[i])
    cmd2 <- paste("bcftools consensus -f $HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
                 "-s", idtab$idchar[i],
                 "largedata/gatk_vcf/JRI20_filtered_snps_annot.bcf.gz",
                 paste0("-o largedata/wgbs_pgenome/", idtab$idchar[i], "/", idtab$idchar[i], ".fa"))
    cat(c(cmd1, cmd2), file=shid, sep="\n", append=FALSE)
}

set_array_job(shid="slurm-script/run_pg.sh", shcode="sh slurm-script/run_pg$SLURM_ARRAY_TASK_ID.sh",
              arrayjobs="1-20", wd=NULL, jobid="pgjob", email="yangjl0930@gmail.com",
              run = c(TRUE, "bigmemm", "2"))

#####################################################################################################
# bcftools view -r 1:1-1000 JRI20_filtered_snps_annot.bcf.gz 
### checking results
"samtools faidx JRIAL2A.fa 1:1-1000"
# Note: heter=> change to alt, multi=>change to dominant alt, missing=> not change

########## preparing genome
for(i in 1:nrow(idtab)){
    shid <- paste0("slurm-script/run_pg", i, ".sh")
    cmd1 <- paste0("bismark_genome_preparation --bowtie2 largedata/wgbs_pgenome/", idtab$idchar[i])
    cat(cmd1, file=shid, sep="\n", append=FALSE)
}    

cmd <- c("module load bismark/0.14.3", "module load bowtie2/2.2.5", "sh slurm-script/run_pg$SLURM_ARRAY_TASK_ID.sh")
set_array_job(shid="slurm-script/run_pg.sh", shcode=cmd,
              arrayjobs="1-20", wd=NULL, jobid="pgjob", email="yangjl0930@gmail.com",
              run = c(TRUE, "bigmemm", "4"))

########### alignment
fq1 <- list.files(path="largedata/wgbs_fq", pattern="R1.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="largedata/wgbs_fq", pattern="R2.fastq.gz$", full.names = TRUE)

#bamfiles <- list.files(path="/group/jrigrp4/BS_teo20/WGBS/BSM", pattern="bam$", full.names = TRUE)
### note: for alignment, "bam" col should not be present.
inputdf <- data.frame(fq1 = fq1,  fq2 = fq2, outbase = gsub(".*/|_.*", "", fq1))
inputdf <- merge(inputdf, idtab[, c("idchar", "pid")], by.x="outbase", by.y="pid")
inputdf$genome <- paste0("$HOME/Documents/Github/methylation/largedata/wgbs_pgenome/", inputdf$idchar)

### AGPv2
run_bismark(inputdf, genome = NULL,
            outdir = "/home/jolyang/Documents/Github/methylation/largedata/wgbs_bismark", 
            N = 1, align = TRUE,
            email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 16))

