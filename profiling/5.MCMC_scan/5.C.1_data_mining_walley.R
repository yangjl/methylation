### Jinliang Yang
### 01-25-2017
### Mining Walley et al., data

# http://science.sciencemag.org/content/353/6301/814/tab-figures-data
info <- read.csv("largedata/Walley_etal/info_maize_2016.csv")

fpkm <- read.csv("largedata/Walley_etal/aag1125_S1_mRNA_FPKM.csv")

fgs <- subset(fpkm, Gene_Set %in% "filtered_set")
set <- fgs[, 1:2]
set$m <- apply(fgs[, 6:28], 1, mean)
set$sd <- apply(fgs[, 6:28], 1, sd)

hist(log(set$m))
median(set$m)
hist(log(set$sd))
median(set$sd)
#brep <- read.csv("largedata/Walley_etal/aag1125_S1_mRNA_bioRep.csv")
dim(set)

write.csv(set, "cache/geneset_rnaseq.csv")


#####################
dnsaf <- read.csv("largedata/Walley_etal/aag1125_S2_dNSAF.csv")
set2 <- dnsaf[, c(1, 50)]
nm <- names(dnsaf)
idx <- which(nm == "Germ.Kernels.2.DAI")
set2$m <- apply(dnsaf[,2:idx], 1, mean)
hist(log(set2$m))

set2$sd <- apply(dnsaf[,2:idx], 1, sd)
hist(log(set2$sd))

write.csv(set2, "cache/geneset_protein.csv")


