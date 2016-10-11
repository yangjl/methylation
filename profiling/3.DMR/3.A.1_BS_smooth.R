### Jinliang Yang
### 10-10-2016
### purpose: BS smooth



library(data.table)
library(bsseq)

meth2 <- fread("largedata/wgbs_smoothed/JRA1_pe.cg.rc")


meth3 <- as.data.frame(meth2)
#meth3$V8 <- meth3$V7 - meth3$V6

methf <- subset(meth3, V5 != "." & V2=="+")



#tmp <- dat[dat$strand == "+",]

BS <- BSseq(chr = meth3$V2, pos = meth3$V3,
            M = as.matrix(meth3$V5, ncol = 1),
            Cov = as.matrix(meth3$V7, ncol = 1), 
            sampleNames = "JRA1_pe_cg")

res <- BSmooth(BS.forward, ns = 70, h = 1000, maxGap = 10^8,
               parallelBy = c("sample", "chromosome"), mc.preschedule = FALSE,
               mc.cores = 1, keep.se = FALSE, verbose = TRUE)

head(getCoverage(res))

head(getMeth(res, type="smooth"))

methf$sm <- getMeth(res, type="smooth")

write.table(methf[, c("chr", "pos", "pos", "sm")], "largedata/chr10.txt", sep="\t", row.names=FALSE, quote=FALSE)

