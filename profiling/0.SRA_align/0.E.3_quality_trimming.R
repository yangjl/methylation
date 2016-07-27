### Jinliang Yang
### check the quality of WGBS data
### 3/22/2016

library("farmeR")

### list files and run QC
fq1 <- list.files(path="largedata/wgbs_fq", pattern="R1.fastq.gz$", full.names=T)
fq2 <- list.files(path="largedata/wgbs_fq", pattern="R2.fastq.gz$", full.names=T)


inputdf <- data.frame(fq1=fq1, fq2=fq2, out1= gsub("_R1", "_R1_trimmed",  fq1), 
                      out2=gsub("_R2", "_R2_trimmed",  fq2))
run_cutadapt(inputdf[1, ], ad1="AGATCGGAAGAGC", ad2="AGATCGGAAGAGC", q=20,
             email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 1))

