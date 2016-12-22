### Jinliang Yang
### 12-06-2016
### determine the cutoff after fitting three normal models

comet <- read.csv("largedata/COMET_CHG/CHG_COMET/JRB2_COMET.csv")
# https://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
library(mixtools)
res <- normalmixEM(x=comet$meth, mu=c(0.01, 0.5, 0.8), k =3)
#lambda The final mixing proportions.
#mu The final mean parameters.
#sigma The final standard deviations. If arbmean = FALSE, then only the smallest standard
#deviation is returned. See scale below

lambda <- summary(res)
mycutoff <- quantile(comet$meth, cumsum(res$lambda))


pdf("graphs/SFig_cutoff_chg_JRB2.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(res, which=1, lwd1=3, col1 ="bisque4")
plot(res, which=2, breaks=50, xlab2="Methylation ratio of COMETs")
lines(density(comet$meth), lty=2, lwd=1.5)
dev.off()


#################
library(mixtools)
library(data.table)

files <- list.files(path="largedata/COMET_CHG/CHG_COMET/", pattern="COMET.csv", full.names = TRUE)
dt <- fread(files[1], data.table=FALSE)
res <- normalmixEM(x=dt$meth, mu=c(0.01, 0.5, 0.8), k =3)
mycutoff <- quantile(dt$meth, cumsum(res$lambda))
out <- data.frame(id=files[1], q1=mycutoff[1], q2=mycutoff[2], q3=mycutoff[3])
for(i in 2:20){
    dt <- fread(files[i], data.table=FALSE)
    res <- normalmixEM(x=dt$meth, mu=c(0.01, 0.5, 0.8), k =3)
    tem <- quantile(dt$meth, cumsum(res$lambda))
    print(tem)
    out2 <- data.frame(id=files[i], q1=tem[1], q2=tem[2], q3=tem[3])
    out <- rbind(out, out2)
}

write.table(out, "cache/chg_cutoff.csv", sep=",", row.names=FALSE, quote=FALSE)

