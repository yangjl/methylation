### Jinliang Yang
### August 4th, 2016

library(farmeR)

files <- list.files(path = "largedata/wgbs_align", pattern = "PE_report.txt", full.names = TRUE)
features <- c("Sequence pairs analysed in total:\t",
              "Number of paired-end alignments with a unique best hit:\t",
              "Sequence pairs with no alignments under any condition:\t",
              "Sequence pairs did not map uniquely:\t",
              #"Sequence pairs which were discarded because genomic sequence could not be extracted:",
              "Total number of C's analysed:\t",
              "Total methylated C's in CpG context:\t",
              "Total methylated C's in CHG context:\t",
              "Total methylated C's in CHH context:\t",
              "Total unmethylated C's in CpG context:\t",
              "Total unmethylated C's in CHG context:\t",
              "Total unmethylated C's in CHH context:\t"
)
res <- get_file2tab(files, features, replace=T )

names(res) <- c("totpe", "hit1", "hit0", "hitN", "totC","mCG", "mCHG", "mCHH", "unCG", "unCHG", "unCHH")
row.names(res) <- gsub(".*/|_PE_.*", "", files) 

res <- as.data.frame(apply(res, 2, as.numeric))
res$seqdp <- (res$totpe)*200/(2500*1e6)
res$mr1 <- with(res, round(hit1/totpe,3)) #uniquely mapped rate
res$mr0 <- with(res, round(hit0/totpe,3)) #non-mapping rate
res$mrN <- with(res, round(hitN/totpe,3)) #multiple-mapping rate

res$cg <- with(res, round(mCG/(mCG+unCG), 3))
res$chg <- with(res, round(mCHG/(mCHG+unCHG), 3))
res$chh <- with(res, round(mCHH/(mCG+unCHH), 3))
res$seqid <- gsub(".*/|_PE_.*", "", files) 

write.table(res, "cache/pg20_bismap_stat.csv", sep=",", row.names=FALSE, quote=FALSE)

######### Plot

pg20 <- read.csv("cache/pg20_bismap_stat.csv")

library(ggplot2)
library(tidyr)

mr <- gather(pg20[, c("seqid", "mr1","mr0", "mrN")], type, reads, 2:4)
mr$type <- factor(mr$type, levels = c("mr1", "mrN", "mr0"))
ratio <- gather(pg20[, c("seqid", "cg","chg", "chh")], type, cr, 2:4)
ratio$type <- factor(ratio$type, levels = c("cg", "chg", "chh"))

###### plot the mapping rate
p1 <- ggplot(mr, aes(x=type, y=reads, fill=type)) +
    geom_violin() +
    theme_bw() +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size=18, face="bold")) +
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0")) +
    scale_x_discrete(labels=c("Unqiue", "Multiple", "Non")) +
    ggtitle("Mapping Rates (N=20)") + xlab("") + 
    ylab("Mapping Rate") + 
    guides(fill=FALSE)
#guides(colour=FALSE, linetype=FALSE)
###### plot the mapping rate
p2 <- ggplot(ratio, aes(x=type, y=cr, fill=type)) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size=18, face="bold")) +
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0")) +
    scale_x_discrete(labels=c("CG", "CHG", "CHH")) +
    ggtitle("Ratio of unmethylated C (N=20)") + xlab("") + 
    ylab("Methylation Ratio") + 
    guides(fill=FALSE)
#guides(colour=FALSE, linetype=FALSE)

source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")
multiplot(p1, p2, cols=2)
