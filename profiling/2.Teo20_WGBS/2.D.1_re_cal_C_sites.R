### Jinliang Yang
### 09-13-2016
### purpose: re-calibrate the C/T and G/A heterozygote sites

library(data.table)

test <- fread("largedata/wgbs_smoothed/JRA1_pe.cg")
chr <- test[V1 == 10]
chr[, snpid := paste(V1, V2, sep="_")]

h <- read.table("largedata/gatk_vcf/JRI20_bi_snps_annot.header", header=T)
snpdt <- fread("largedata/gatk_vcf/JRI20_bi_snps_annot.txt")
names(snpdt) <- names(h)

sid <- "JRIAL2A"
sub <- snpdt[, c("chr", "pos", "ref", "alt", sid),  with=FALSE]
sub <- sub[chr == 10]
names(sub)[5] <- "SAMPLE"
sub <- sub[SAMPLE == "Y" | SAMPLE == "R"]
#sub$snpid <- paste(sub$chr, sub$pos, sep="_")

sub[, snpid := paste(chr, pos, sep="_")]

out <- merge(chr, sub, by.x="snpid", by.y="snpid", all.x=TRUE)
out[, tot := V4 + V5]
out[!is.na(SAMPLE), tot := ceiling((V4 + V5)/2)]
out[V4 > tot & SAMPLE == "Y", tot := V4]

