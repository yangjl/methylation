### Jinliang Yang
### 03-20-2017
### exploring the function of gbm

library("data.table")

res <- read.csv("cache/stat_exon_mean_var.csv")


tpm <- fread("largedata/tpm_9runs.csv", data.table=FALSE)
tpm$geneid <- gsub(".exon.*", "", tpm$target_id)
tpm$geneid <- gsub("_E.*", "", tpm$geneid)


d <- merge(res, tpm, by="geneid")
d$exon <- gsub(".*exon|.*E", "", d$target_id)
d$exon <- as.numeric(as.character(d$exon))




library("beanplot")

myd <- d[, c("geneid", "mm", "exon", "SRR957418")]
names(myd)[4] <- "tpm"
    


par(mfrow=c(1,2))
 
d0 <- subset(myd, mm < 0.1 & exon < 10)
out0 <- subset(d0, tpm > 2)
dim(out0)
beanplot(log(tpm) ~ exon, data = out0, kernel="cosine", ll = 0.04, cex=fs, side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"))
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)

d1 <- subset(myd, mm >= 0.5 & exon < 10)
out1 <- subset(d1, tpm > 2)
dim(out1)
beanplot(log(tpm) ~ exon, data = out1, kernel="cosine", ll = 0.04, cex=fs, side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"))
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)








par(mfrow=c(2,2))
d1 <- subset(d, mm > 0.6 & exon < 5)
out1 <- subset(d1, SRR957418 > 1)
dim(out1)
beanplot(log2(SRR957418) ~ exon, data = out1, kernel="cosine", ll = 0.04, cex=fs, side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"))
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)


tpm <- subset(tpm, geneid %in% gbM$geneid)



####
gbM <- subset(res, mm > 0.6)

gff$geneid <- gsub(";.*|_.*", "", gff$attribute)
gff$geneid <- gsub(".*=", "", gff$geneid)
gff <- subset(gff, geneid %in% gbM$geneid)






