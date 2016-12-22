### Jinliang Yang
### 12-06-2016
### determine the cutoff after fitting three normal models

comet <- read.csv("largedata/COMET/CG_COMET/JRB1_CG_COMET.csv")

# https://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
library(mixtools)
res <- normalmixEM(x=comet$meth, mu=c(0.01, 0.5, 0.8), k =3)
#lambda The final mixing proportions.
#mu The final mean parameters.
#sigma The final standard deviations. If arbmean = FALSE, then only the smallest standard
#deviation is returned. See scale below

lambda <- summary(res)

pdf("graphs/SFig_cutoff.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(res, which=1, lwd1=3, col1 ="bisque4")
plot(res, which=2, breaks=50, xlab2="Methylation ratio of COMETs")
lines(density(comet$meth), lty=2, lwd=1.5)
dev.off()


#################
library(mixtools)

files <- list.files(path="largedata/COMET_CHG/CHG_COMET/", pattern="COMET.csv", full.names = TRUE)
dt <- fread(files[i], data.table=FALSE)
res <- normalmixEM(x=dt$meth, mu=c(0.01, 0.5, 0.8), k =3)
mycutoff <- quantile(dt$meth, cumsum(res$lambda))

for(i in 12:20){
    dt <- fread(files[i], data.table=FALSE)
    res <- normalmixEM(x=dt$meth, mu=c(0.01, 0.5, 0.8), k =3)
    tem <- quantile(dt$meth, cumsum(res$lambda))
    mycutoff <- rbind(mycutoff, tem)
}

write.table(mycutoff, "cache/chg_cutoff.csv", sep=",", row.names=FALSE, quote=FALSE)

