### Jinliang Yang
### 03-29-2017
### B73 five tissues

#########
library("GenomicRanges")
library("data.table")
library(plyr)
source("lib/find_overlap_win_gff.R")

get_gbm <- function(){
    
    ####### AGPv4 GFF, but use only 1-to-1 from v3 vs. v4
    gff <- read.delim("~/dbcenter/AGP/AGPv4/Zea_mays.AGPv4.34.gff3", comment.char = "#", header=FALSE)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    ### read v3=>v4 genes
    gene34 <- read.table("~/dbcenter/AGP/AGPv4/maize.v3TOv4.geneIDhistory.txt", header=TRUE)
    gene11 <- subset(gene34, type %in% "1-to-1")
    
    
    exon <- subset(gff, feature %in% "exon")
    exon$txid <- gsub(".*transcript:|;Name=.*", "", exon$attribute)
    idx <- grep("T001$|T01$", exon$txid)
    
    exon1 <- exon[idx, ]
    exon1$exonid <- gsub(".*exon", "exon", exon1$attribute)
    exon1$exonid <- gsub(";.*", "", exon1$exonid)
    exon1$geneid <- gsub("_.*", "", exon1$txid)
    
    exon1 <- subset(exon1, geneid %in% gene11$V4_gene)
    length(unique(exon1$geneid))
    # 28390
    
    ### readin data of B73 with five tissues
    b73 <- fread("largedata/ZmV4_w100.fivetissues.header.bed", header=FALSE, skip=1L, data.table=FALSE)
    h <- read.delim("largedata/ZmV4_w100.fivetissues.header.bed", header=TRUE, nrows=3)
    h <- h[, -4]
    names(b73) <- names(h)
    
    idx <- grep("Cratio", names(b73))
    for(i in idx){
        
        cr <- names(b73)[i]
        df <- b73[, c("chr", "start", "end", cr)]
        names(df)[4] <- "ratio"
        
        message(sprintf("###>>> working on [ %s ] ...", cr))
        m <- find_overlap_win_gff(df, gff=exon1)
        mtx <- m[[1]]
        mexon <- m[[2]]
        
        write.table(mtx, paste0("largedata/B73_meth/body_", cr,".csv"), sep=",", row.names=FALSE, quote=FALSE)
        write.table(mexon, paste0("largedata/B73_meth/exon_", cr, ".csv"), sep=",", row.names=FALSE, quote=FALSE)
    }
}

#####
get_gbm()
###>>> working on [ Cratio_B73_all3_CG ] ...
###>>> working on [ Cratio_B73_all3_CHG ] ...
###>>> working on [ Cratio_B73_all3_CHH ] ...
###>>> working on [ Cratio_B73_anther_CG ] ...
###>>> working on [ Cratio_B73_anther_CHG ] ...
###>>> working on [ Cratio_B73_anther_CHH ] ...
###>>> working on [ Cratio_B73_earshoot_CG ] ...
###>>> working on [ Cratio_B73_earshoot_CHG ] ...
###>>> working on [ Cratio_B73_earshoot_CHH ] ...
###>>> working on [ Cratio_B73_flag_leaf_CG ] ...
###>>> working on [ Cratio_B73_flag_leaf_CHG ] ...
###>>> working on [ Cratio_B73_flag_leaf_CHH ] ...
###>>> working on [ Cratio_B73_SAM_CG ] ...
###>>> working on [ Cratio_B73_SAM_CHG ] ...
###>>> working on [ Cratio_B73_SAM_CHH ] ...


collect_body <- function(pwd, pt="body", context="CG"){
    bfiles <- list.files(path=pwd, pattern=pt, full.names = TRUE)
    
    cgidx <- grep(context, bfiles)
    out <- read.csv(bfiles[cgidx[1]])
    out <- out[, 1:2]
    names(out)[2] <- gsub(".*\\/|.csv", "", bfiles[cgidx[1]])
    
    for(i in 2:length(cgidx)){
        tem <- read.csv(bfiles[cgidx[i]])
        tem <- tem[, 1:2]
        names(tem)[2] <- gsub(".*\\/|.csv", "", bfiles[cgidx[i]])
        out <- merge(out, tem, by="txid")
    }
    out$mean <- apply(out[, -1], 1, function(x) mean(x, na.rm=TRUE))
    out$var <- apply(out[, -1], 1, function(x) var(x, na.rm=TRUE))
    return(out)
}

########
out1 <- collect_body(pwd="largedata/B73_meth", pt="body", context="CG")
out1 <- out1[order(out1$var, decreasing=T), ]
write.table(out1, "largedata/B73_5tissues_body_cg.csv", sep=",", row.names=FALSE, quote=FALSE)

out2 <- collect_body(pwd="largedata/B73_meth", pt="body", context="CHG")
out2 <- out2[order(out2$var, decreasing=T), ]
write.table(out2, "largedata/B73_5tissues_body_chg.csv", sep=",", row.names=FALSE, quote=FALSE)

out3 <- collect_body(pwd="largedata/B73_meth", pt="body", context="CHH")
out3 <- out3[order(out3$var, decreasing=T), ]
write.table(out3, "largedata/B73_5tissues_body_chh.csv", sep=",", row.names=FALSE, quote=FALSE)

g1 <- as.character(subset(out1, var > 0.01)$txid)
g2 <- as.character(subset(out2, var > 0.01)$txid)
g3 <- as.character(subset(out3, var > 0.001)$txid)
length(g1) #365
length(g2) #227
length(g3) #228
sum(g1 %in% g2) #79
sum(g1 %in% g3) #40
sum(g2 %in% g3) #45

sum(g1[(g1 %in% g2)] %in% g3) #23

g1[g1[(g1 %in% g2)] %in% g3]


collect_exon <- function(pwd, pt="body", context="CG"){
    bfiles <- list.files(path=pwd, pattern=pt, full.names = TRUE)
    
    idx <- grep(context, bfiles)
    out <- read.csv(bfiles[idx[1]])
    out <- out[, 1:2]
    names(out)[2] <- gsub(".*\\/|.csv", "", bfiles[idx[1]])
    
    for(i in 2:length(cgidx)){
        tem <- read.csv(bfiles[cgidx[i]])
        tem <- tem[, 1:2]
        names(tem)[2] <- gsub(".*\\/|.csv", "", bfiles[cgidx[i]])
        out <- merge(out, tem, by="exonid")
    }
    out$mean <- apply(out[, -1], 1, function(x) mean(x, na.rm=TRUE))
    out$var <- apply(out[, -1], 1, function(x) var(x, na.rm=TRUE))
    return(out)
}
########
out1 <- collect_exon(pwd="largedata/B73_meth", pt="exon", context="CG")
out1 <- out1[order(out1$var, decreasing=T), ]
write.table(out1, "largedata/B73_5tissues_exon_cg.csv", sep=",", row.names=FALSE, quote=FALSE)

out2 <- collect_exon(pwd="largedata/B73_meth", pt="exon", context="CHG")
out2 <- out2[order(out2$var, decreasing=T), ]
write.table(out2, "largedata/B73_5tissues_exon_chg.csv", sep=",", row.names=FALSE, quote=FALSE)

out3 <- collect_exon(pwd="largedata/B73_meth", pt="exon", context="CHH")
out3 <- out3[order(out3$var, decreasing=T), ]
write.table(out3, "largedata/B73_5tissues_exon_chh.csv", sep=",", row.names=FALSE, quote=FALSE)

out1$exon <- as.character(gsub(".*exon", "", out1$exonid))
table(out1$exon)

out1 <- read.csv("largedata/B73_5tissues_exon_cg.csv")
out2 <- read.csv("largedata/B73_5tissues_exon_chg.csv")
out3 <- read.csv("largedata/B73_5tissues_exon_chh.csv")

g1 <- as.character(subset(out1, var > 0.1)$exonid)
g2 <- as.character(subset(out2, var > 0.15)$exonid)
g3 <- as.character(subset(out3, var > 0.15)$exonid)
length(g1) #1203
length(g2) #1793
length(g3) #2574
sum(g1 %in% g2) #155
sum(g1 %in% g3) #160
sum(g2 %in% g3) #1449

g1[(g1 %in% g2)]
g2[(g1 %in% g2)]
sum(g1[(g1 %in% g2)] %in% g3) #23

op <- data.frame(exonid=g1[(g1 %in% g2)], exon=1)
op$exon <- as.character(gsub(".*exon", "", op$exonid))

par(mfrow=c(1,2))
plot(table(out1$exon)/nrow(out1))
plot(table(op$exon)/nrow(op))

length(g1[g1[(g1 %in% g2)] %in% g3]/nrow)
