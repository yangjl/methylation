

library("data.table")

infile="largedata/wgbs_smoothed/JRA1_pe.cg.rc"

meth <- fread(infile, data.table=FALSE)

meth <- subset(meth, V2 == 10 & V4 == "+")
meth$ratio <- round(meth$V5/meth$V7, 4)

write.table(meth[, c("V3", "ratio")], "largedata/sub_cg_chr10.csv", sep=",", row.names=FALSE, quote=FALSE )

sub <- read.csv("largedata/sub_cg_chr10.csv")


names(sub)[1] <- "pos1"
sub <- sub[order(sub$pos1),]
sub <- sub[sub$ratio > 0.8, ]
sub$pos2 <- c(sub$pos1[-1], sub$pos1[nrow(sub)])
sub$dis <- sub$pos2 - sub$pos1
mydis <- sub$dis
mydis <- mydis[mydis < 100]
hist(mydis)




###>>> CG enriched genes
#### plotting
library(data.table)
out <- fread("largedata/kmer_count.csv", sep=",", header=FALSE)
h <- read.csv("largedata/kmer_count.csv", sep=",", nrow=5, header=TRUE)
out[, 1:10, with=FALSE]
out <- as.data.frame(out)
names(out) <- c("geneid", names(h))

out2 <- out[, c("geneid", "A", "T", "C", "G")]
out2$cgfrq <- (out2$C + out2$G)/(out2$A + out2$T + out2$C + out2$G)
out2 <- out2[order(out2$cgfrq, decreasing = T), ]
idx <- grep("^10.", as.character(out2$geneid))
head(out2[idx, ])

out3 <- out2[idx,]
out3$geneid <- gsub("^10.", "", out3$geneid)

gene <- get_feature(gff="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", features="gene")
subgene <- merge(out3, gene, by.x="geneid", by.y="attribute")

subgene <- subgene[order(subgene$cgfrq, decreasing = T), ]
head(subgene[, -2:-5], 30)

###>>>>
reg <- c(143025209, 143085610)
myd <- subset(sub, pos1 > reg[1] & pos2 < reg[2] & ratio > 0.8)

plot(x=myd$pos1, y=myd$ratio, type="h")


out <- myd$dis   
out <- out[out < 500]    
table(out)    
hist(table(out), type="h")

