### Jinliang Yang
### 10-16-2016
### purpose: get SFS for features


get_feature_sfs <- function(cometfile="largedata/COMET/CG_COMET/chr10_comet_blocks_0.33.csv", 
                            feature,
                          outfile="cache/sfs_cg_comet.csv",
                          cols=3:22){
    comet <- fread(cometfile, data.table=FALSE)
    
    ### removing missing data
    comet$miss <- apply(comet[, cols], 1, function(x){
        return(sum(x < 0))
    })
    comet <- subset(comet, miss == 0)
    
    ###
    ###########
    chr <- comet[, 1:2]
    chr$start <- as.numeric(as.character(gsub("_.*", "", comet$bid)))
    chr$end <- as.numeric(as.character(gsub(".*_", "", comet$bid)))
    gr1 = with(chr, GRanges(seqnames=chr, IRanges(start=start, end=end), strand="+", bid=bid))
    
    feature <- gr0
    ex1 = findOverlaps(feature, gr1)
    out <- gr1[subjectHits(ex1)]
    
    mycomet <- subset(comet, bid %in% mcols(out)$bid)
    
    ###########
    f0 <- apply(mycomet[, cols], 1, function(x){
        n2 <- sum(x == 2)
        n1 <- sum(x == 1)
        return(2*n2 + n1)
    })
    sfs <- table(f0)
    write.table(sfs, outfile, sep=",", row.names=FALSE, quote=FALSE)
    return(sfs)
}



library("GenomicRanges")
library("data.table")

ob <- load("largedata/AGPv2_features.RData")


#### gene
f[['gene']]$seqname <- paste0("chr", f[['gene']]$seqname)
gr0 = with(f[['gene']], GRanges(seqnames=seqname, IRanges(start=start, end=end), strand=strand, geneid=attribute ))

res1 <- get_feature_sfs(cometfile="largedata/COMET/CG_COMET/chr10_comet_blocks_0.33.csv", 
                feature=gr0, outfile="cache/sfs_cg_comet_0.33_gene.csv", cols=3:22)

#### exon
f[['exon']]$seqname <- paste0("chr", f[['exon']]$seqname)
gr0 = with(f[['exon']], GRanges(seqnames=seqname, IRanges(start=start, end=end), strand=strand, geneid=attribute ))

res2 <- get_feature_sfs(cometfile="largedata/COMET/CG_COMET/chr10_comet_blocks_0.33.csv", 
                       feature=gr0, outfile="cache/sfs_cg_comet_0.33_exon.csv", cols=3:22)
#### intron
f[['intron']]$seqname <- paste0("chr", f[['intron']]$seqname)
gr0 = with(f[['intron']], GRanges(seqnames=seqname, IRanges(start=start, end=end), strand=strand, geneid=attribute ))

res3 <- get_feature_sfs(cometfile="largedata/COMET/CG_COMET/chr10_comet_blocks_0.33.csv", 
                       feature=gr0, outfile="cache/sfs_cg_comet_0.33_intron.csv", cols=3:22)
