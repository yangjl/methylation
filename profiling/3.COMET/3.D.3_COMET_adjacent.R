### Jinliang Yang
### 11-30-2016
### study the feature associated with different methylation classes.

###########
library(data.table)
library(GenomicRanges)

get_overlap <- function(myquery=vout0, gff, repeats){
    
    ### gene features: exon
    exon <- subset(gff, feature %in% "exon")
    grexon <- with(exon, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    
    ### intron
    intron <- subset(gff, feature %in% "intron")
    grin <- with(intron, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    
    ### upstream 1kb regions
    p <- subset(gff, feature %in% "gene" & strand %in% "+")
    p$end <- p$start - 1
    p$start <- p$start - 1001
    
    m <- subset(gff, feature %in% "gene" & strand %in% "-")
    m$start <- m$end + 1
    m$end <- m$end + 1001
    
    onek <- rbind(p, m)
    gronek <- with(onek, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    
    ### Class I
    c1 <- subset(repeats, class %in% "Class I Retroelements")
    grc1 <- with(c1, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    
    ### Class II/III
    c2 <- subset(repeats, class %in% "Class II/III Transposable Ele")
    grc2 <- with(c2, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    
    ### Other TE
    c3 <- subset(repeats, class %in% "Other")
    grc3 <- with(c3, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    
    ###################################################
    myquery$chr <- gsub("chr", "", myquery$chr)
    
    myq <- with(myquery, GRanges(seqnames=chr, IRanges(start=start, end=end)))
    
    ###which features from the ‘query’ overlap which features in the subject’
    #ex1 <- findOverlaps(query=gr1, subject=mysub)
    #ranges(gr0)[queryHits(ex1)] = ranges(gr1)[subjectHits(ex1)]
    
    #mysub$strand <- "*"
    int1 <- intersect(myq, grexon)
    len1 <- as.data.frame(int1)
    
    int2 <- intersect(myq, grin)
    len2 <- as.data.frame(int2)
    
    int3 <- intersect(myq, gronek)
    len3 <- as.data.frame(int3)
    
    int4 <- intersect(myq, grc1)
    len4 <- as.data.frame(int4)
    
    int5 <- intersect(myq, grc2)
    len5 <- as.data.frame(int5)
    
    int6 <- intersect(myq, grc3)
    len6 <- as.data.frame(int6)
    
    
    res <- data.frame(exon=sum(len1$width),
                      intron=sum(len2$width),
                      onek=sum(len3$width),
                      c1=sum(len4$width),
                      c2=sum(len5$width),
                      c3=sum(len6$width)  )
    
    return(res)
}

######
ob <- load("largedata/chr10_comet.RData")
# vout0:low; vout2:high; vout3:variable


gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")

repeats <- fread("~/dbcenter/AGP/AGPv2/repeats/ZmB73_5a_MIPS_repeats.gff", header=FALSE, data.table=FALSE)
names(repeats) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
repeats$class <- gsub(".*type=|;name=.*", "", repeats$attribute)

############
res0 <- get_overlap(myquery=vout0, gff, repeats)
res0$type <- "low"
res2 <- get_overlap(myquery=vout2, gff, repeats)
res2$type <- "high"
res3 <- get_overlap(myquery=vout3, gff, repeats)
res3$type <- "variable"
    
res <- rbind(res0, res2, res3)
res$tot <- apply(res[, 1:6], 1, sum)
write.table(res, "cache/feature_chr10.csv", sep=",", row.names=FALSE, quote=FALSE)

out <- res
out$exon <- with(out, exon/tot)
out$intron <- with(out, intron/tot)
out$onek <- with(out, onek/tot)
out$c1 <- with(out, c1/tot)
out$c2 <- with(out, c2/tot)
out$c3 <- with(out, c3/tot)




  