### Jinliang Yang
### August 4th, 2016

library(farmeR)

files <- list.files(path = "largedata/bismark", pattern = "PE_report.txt", full.names = TRUE)
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

write.table(res, "cache/comp3_bismap_stat.csv", sep=",", row.names=FALSE, quote=FALSE)

######### Plot

comp <- read.csv("cache/comp3_bismap_stat.csv")

library(ggplot2)
library(tidyr)

lres <- gather(comp[, c("seqid", "mr1","mr0", "mrN")], type, reads, 2:4)


#lres <- lres[order(lres$type, lres$type, decreasing = TRUE),]
lres$type <- factor(lres$type, levels = c("mr1", "mrN", "mr0"), 
                    labels=c("Unqiue", "Multiple", "Non"), ordered=TRUE)

theme_set(theme_grey(base_size = 18)) 
s <- ggplot(lres, aes(x=seqid, y=reads, fill = type)) + 
    #opts(axis.text.x=theme_text(angle=90)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="", y="Mapping Rate", fill="Type") +
    scale_x_discrete(labels=c("B73", "B73 -> N", "B73 -> pseduoRef"))
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12)
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))
s

