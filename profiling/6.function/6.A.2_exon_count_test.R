### Jinliang Yang
### 03-20-2017
### exploring the function of gbm

library("data.table")
library(tidyr)

res <- read.csv("cache/stat_exon_mean_var.csv")


tpm <- fread("largedata/tpm_9runs.csv", data.table=FALSE)
tpm$geneid <- gsub(".exon.*", "", tpm$target_id)
tpm$geneid <- gsub("_E.*", "", tpm$geneid)


d <- merge(res, tpm, by="geneid")
d$exon <- gsub(".*exon|.*E", "", d$target_id)
d$exon <- as.numeric(as.character(d$exon))


out <- gather(d, sample, exp, 8:16)

rtb <- read.delim("data/SraRunTable_walley_etal_2016.txt", header=TRUE)
#c('submission','study','sample','experiment')
sum(rtb$MBases_l)

sra <- data.frame(SRR=rtb$Run_s, sid=rtb$Sample_Name_s, pid=rtb$tissue_s, gid=rtb$genotype_s)

out2 <- merge(out, sra[, c("SRR", "pid")], by.x="sample", by.y="SRR")
out2$uid <- paste(out2$sample, out2$geneid, sep="_")

exon1 <- subset(out2, exon == 1)

exon2 <- subset(out2, exon == 2)


exon1_2 <- merge(exon1, exon2[, c("uid", "exp")], by="uid")


exon1_2 <- subset(exon1_2, exp.x > 1 & exp.y > 1)

exon1_2$exp21 <- exon1_2$exp.y - exon1_2$exp.x
exon1_2 <- subset(exon1_2, exp21 > 0)

plot(exp21 ~ mm, data=exon1_2, cex=0.3, pch=16)



p1 <- subset(exon1_2, exp21 > 0)
par(mfrow=c(1,2))
plot(log(exp21) ~ mm, data=p1, cex=0.3, pch=16)
abline(h = 0)

p2 <- subset(exon1_2, exp21 < 0)
plot(log(abs(exp21)) ~ mm, data=p2, cex=0.3, pch=16)
abline(h = 0)


tem <- subset(exon1_2, log(exp21) < -10)

ngbm <- subset(exon1_2, mm < 0.1)
hist(ngbm$exp.x - ngbm$exp.y, breaks=100)
hist()




########
exon1_2 <- subset(exon1_2, exp.y > 1 & exp.x > 1)
exon1_2$exp21 <- log(exon1_2$exp.y) - log(exon1_2$exp.x)
#exon1_2 <- subset(exon1_2, abs(exp12) > 0)
fit <- lm(exp21 ~ mm + sample + pid, data=exon1_2)

fit <- glm( mm ~ exp21 + sample + pid, data=exon1_2, family="binomial")
plot(fit)

exon1_2 <- subset(exon1_2, abs(exp12) > 0)
plot(exp21 ~ mm, data=exon1_2, cex=0.2, pch=16)


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
beanplot(log(SRR957418) ~ exon, data = out1, kernel="cosine", ll = 0.04, cex=fs, side = "no", cut=10,
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"))
#axis(side =1, at =1:4, labels =c("25", "50", "75", "100"), cex.axis=fs)


tpm <- subset(tpm, geneid %in% gbM$geneid)



####
gbM <- subset(res, mm > 0.6)

gff$geneid <- gsub(";.*|_.*", "", gff$attribute)
gff$geneid <- gsub(".*=", "", gff$geneid)
gff <- subset(gff, geneid %in% gbM$geneid)






