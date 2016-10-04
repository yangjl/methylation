### Jinliang Yang
### 10-04-2016
### purpose: BS smooth



library(data.table)
library(bsseq)

meth2 <- fread("largedata/wgbs_smoothed/JRA1_pe.cg.relc")


meth3 <- as.data.frame(meth2)

methf <- subset(meth3, V5 != "." & V2=="+")



#tmp <- dat[dat$strand == "+",]

BS.forward <- BSseq(chr = meth3$V2, pos = meth3$V3,
                    M = as.matrix(methf$V4, ncol = 1),
                    Cov = as.matrix(methf$V5, ncol = 1), sampleNames = "forward")

res <- BSmooth(BS.forward, ns = 70, h = 1000, maxGap = 10^8,
               parallelBy = c("sample", "chromosome"), mc.preschedule = FALSE,
               mc.cores = 1, keep.se = FALSE, verbose = TRUE)

head(getCoverage(res))

head(getMeth(res, type="smooth"))

methf$sm <- getMeth(res, type="smooth")

write.table(methf[, c("chr", "pos", "pos", "sm")], "largedata/chr10.txt", sep="\t", row.names=FALSE, quote=FALSE)

