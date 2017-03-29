### Jinliang Yang
### 03-29-2017
### gene body ratio of B73




mtx <- read.csv("largedata/B73_meth/body_Cratio_B73_all3_CG.csv")
mtx <- subset(mtx, !is.na(mean))
hist(mtx$mean, breaks=100)

nrow(subset(mtx, mean > 0.6))


mexon <- read.csv("largedata/B73_meth/exon_Cratio_B73_all3_CG.csv")
mexon <- subset(mexon, !is.na(mean))
mexon$exonid <- gsub(".*exon", "", mexon$exonid)

library("beanplot")
par(mfrow=c(1,2))


mexon6 <- subset(mexon, V4_gene %in% subset(mtx, mean > -1)$V4_gene)
mexon6$exonid <- as.numeric(as.character(mexon6$exonid))
mexon6 <- subset(mexon6,  exonid < 11)

beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"),
         xlab="Exon#", ylab="CG methylation", main="B73 CG methylation")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)


mexon6 <- subset(mexon, V4_gene %in% subset(mtx, mean > 0.8)$V4_gene)
mexon6$exonid <- as.numeric(as.character(mexon6$exonid))
mexon6 <- subset(mexon6,  exonid < 11)

beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"),
         xlab="Exon#", ylab="CG methylation", main="Gene body methylated Genes")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)

mexon6 <- subset(mexon, V4_gene %in% subset(mtx, mean < 0.8)$V4_gene)
mexon6$exonid <- as.numeric(as.character(mexon6$exonid))
mexon6 <- subset(mexon6,  exonid < 11)

beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"),
         xlab="Exon#", ylab="CG methylation", main="Non-gbM Genes")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)



mtx <- read.csv("largedata/B73_meth/body_Cratio_B73_all3_CHG.csv")
mtx$geneid <- gsub("_.*", "", mtx$txid)
mtx <- subset(mtx, !is.na(mean))
hist(mtx$mean, breaks=100)

nrow(subset(mtx, mean > 0.6))


mexon <- read.csv("largedata/B73_meth/exon_Cratio_B73_all3_CHG.csv")
mexon$geneid <- gsub("_.*", "", mexon$exonid)
mexon <- subset(mexon, !is.na(mean))
mexon$exonid <- gsub(".*exon", "", mexon$exonid)
mexon$exonid <- as.numeric(as.character(mexon$exonid))




######################
library("beanplot")
#par(mfrow=c(1,2))
# col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"),

mexon6 <- subset(mexon,  exonid < 11)
mexon6 <- subset(mexon6, geneid %in% subset(mtx, mean > -1)$geneid)
beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, 
         xlab="Exon#", ylab="CHG methylation", main="B73 CHG methylation")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)

mexon6 <- subset(mexon,  exonid < 11)
mexon6 <- subset(mexon6, geneid %in% subset(mtx, mean > 0.8)$geneid)
beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("antiquewhite3"),
         xlab="Exon#", ylab="CHG methylation", main="Gene body methylated Genes")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)

mexon6 <- subset(mexon,  exonid < 11)
mexon6 <- subset(mexon6, geneid %in% subset(mtx, mean < 0.8)$geneid)
beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("antiquewhite3"),
         xlab="Exon#", ylab="CHG methylation", main="Non-gbM Genes")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)








mtx <- read.csv("largedata/B73_meth/body_Cratio_B73_all3_CHH.csv")
mtx$geneid <- gsub("_.*", "", mtx$txid)
mtx <- subset(mtx, !is.na(mean))
#hist(mtx$mean, breaks=100)


mexon <- read.csv("largedata/B73_meth/exon_Cratio_B73_all3_CHH.csv")
mexon$geneid <- gsub("_.*", "", mexon$exonid)
mexon <- subset(mexon, !is.na(mean))
mexon$exonid <- gsub(".*exon", "", mexon$exonid)
mexon$exonid <- as.numeric(as.character(mexon$exonid))


######################
library("beanplot")
#par(mfrow=c(1,2))
# col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"),
par(mfrow=c(3,1))

mexon6 <- subset(mexon,  exonid < 11)
mexon6 <- subset(mexon6, geneid %in% subset(mtx, mean > -1)$geneid)
beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, 
         xlab="Exon#", ylab="CHG methylation", main="B73 CHG methylation")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)

mexon6 <- subset(mexon,  exonid < 11)
mexon6 <- subset(mexon6, geneid %in% subset(mtx, mean > 0.1)$geneid)
beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("antiquewhite3"),
         xlab="Exon#", ylab="CHG methylation", main="Gene body methylated Genes")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)

mexon6 <- subset(mexon,  exonid < 11)
mexon6 <- subset(mexon6, geneid %in% subset(mtx, mean < 0.1)$geneid)
beanplot(mean ~ exonid, data = mexon6, kernel="cosine", overalline="median", ll = 0.005,  side = "no", cut=10,
         border = NA, col=list("antiquewhite3"),
         xlab="Exon#", ylab="CHG methylation", main="Non-gbM Genes")
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)