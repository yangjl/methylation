### Jinliang Yang
### 8/11/2016
### Counting kmer frq for each genes

####### get gene id
geneid <- read.table("/home/jolyang/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt", header=T)
geneid <- subset(geneid, is_canonical == "yes")


feature <- geneid[, c("chromosome", "transcript_start", "transcript_end", "gene_id")]
names(feature) <- c("chr", "start", "end", "geneid")
feature$chr <- gsub("chr", "", feature$chr)

feature <- subset(feature, chr %in% 1:10)

library("Biostrings")
library("pseudoRef")
fa <- readDNAStringSet(filepath = "~/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa", format="fasta")

res <- kcount(fa, feature, kmers=1:7)
out <- t(res)
write.table(out, "largedata/kmer_count.csv", sep=",", quote=FALSE)


#### plotting
library(data.table)
out <- fread("largedata/kmer_count.csv", sep=",", header=TRUE)
## get mean count
mc <- apply(out, 2, mean)






