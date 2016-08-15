### Jinliang Yang
### 8/15/2016
### ml

library(data.table)

geneid <- fread("/home/jolyang/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt", header=T, data.table=FALSE)
geneid <- subset(geneid, is_canonical == "yes")

feature <- geneid[, c("chromosome", "transcript_start", "transcript_end", "gene_id")]
names(feature) <- c("chr", "start", "end", "geneid")
feature$chr <- gsub("chr", "", feature$chr)

feature <- subset(feature, chr %in% 1:10)
feature$len <- feature$end - feature$start

###################
gbm <- read.csv("cache/stat_exon_mean_var.csv")
gbm$mb <- 1
gbm[gbm$mm < 0.6, ]$mb <- 0
gbm <- merge(feature[, c("geneid", "len")], gbm, by="geneid")

###################
out <- fread("largedata/kmer_count.csv", sep=",", header=FALSE)
h <- read.csv("largedata/kmer_count.csv", sep=",", nrow=5, header=TRUE)
out[, 1:10, with=FALSE]
out <- as.data.frame(out)
names(out) <- c("geneid", names(h))

out$geneid <- gsub("^.{1,3}\\.", "", out$geneid)


df <- merge(gbm, out, by="geneid")
df[, -1:-5] <- round(df[, -1:-5]/df$len,4)

### clean up data
rm(list=c("out", "feature", "geneid", "h"))

#################################################################################
library(glmnet)
#require(doMC)
#registerDoMC(cores=6)

tb <- data.frame(na=names(df), kmer=1)
tb$kmer <- nchar(as.character(tb$na))

### train binomial model
#mfit <- glmnet(x= as.matrix(df[, idx]), y=df$mb, family = "binomial")
idx <- 5+which(tb$kmer[-1:-5] <= 3)
cvfit = cv.glmnet(x= as.matrix(df[, idx]), y=df$mb, family = "binomial", type.measure = "class", parallel=FALSE)
cvfit$lambda.min
head(sort(coef(cvfit, s = "lambda.min")))

save(file="largedata/cvfit.RData", list="cvfit")
#res <- coef(cvfit, s = "lambda.min")


ob <- load("largedata/cvfit.RData")
res <- as.matrix(coef(cvfit, s = "lambda.min"))

head(res[order(res[,1]),], 30)

print(cvfit)
plot(cvfit)
