### Jinliang Yang
### July 25th, 2016

library(farmeR)

idtab <- read.csv("data/teo20_ids.csv")
idtab$idchar <- gsub("-", "", idtab$idchar)
idtab$pid <- paste0("JR", idtab$plate)
########### alignment
fq1 <- list.files(path="largedata/wgbs_fq", pattern="R1_trimmed.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="largedata/wgbs_fq", pattern="R2_trimmed.fastq.gz$", full.names = TRUE)

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

