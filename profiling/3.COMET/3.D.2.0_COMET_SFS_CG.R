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


############## get annotation
gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")

repeats <- fread("~/dbcenter/AGP/AGPv2/repeats/ZmB73_5a_MIPS_repeats.gff", header=FALSE, data.table=FALSE)
names(repeats) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
repeats$class <- gsub(".*type=|;name=.*", "", repeats$attribute)
################ get SFS of CG


######
library("data.table")
df <- read.csv("largedata/lcache/SFS_comet_blocks_CG.csv")
df$chr <- gsub("chr", "", df$chr)

out <- data.frame()
for(sitei in 0:40){
    mydf <- subset(df, sfs == sitei)
    tmp <- get_overlap(myquery=mydf, gff, repeats)
    tmp$site <- sitei
    out <- rbind(out, tmp)
}

write.table(out, "cache/SFS_comet_features_CG.csv", sep=",", row.names=FALSE, quote=FALSE)



