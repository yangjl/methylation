
### find COMET overlap with features and then run MCMC for the feature SFS
run_overlap_MCMC <- function(df, gff, fea, runid, outdir){
    # df:
    # gff: genomic feature
    # feature: vector [exon, intron, gene, Class I Retroelements, etc.]
    # outdir: "largedata/lcache/"
    
    ###################################################
    ### GR of COMET
    #myquery$chr <- gsub("chr", "", myquery$chr)
    grc <- with(df, GRanges(seqnames=chr, IRanges(start=start, end=end), cid=cid))
    
    if(fea != "gene"){
        ### gene features: exon
        f <- subset(gff, feature %in% fea)
        grf <- with(f, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
        
    }else if (fea == "gene"){
        ### upstream 1kb regions
        p <- subset(gff, feature %in% "gene" & strand %in% "+")
        p$end <- p$start - 1
        p$start <- p$start - 1001
        
        m <- subset(gff, feature %in% "gene" & strand %in% "-")
        m$start <- m$end + 1
        m$end <- m$end + 1001
        
        onek <- rbind(p, m)
        grf <- with(onek, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
    }
    
    out <- get_overlap_sfs(grf, grc, df)
    # If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
    res <- MCMCBC(my_sfs=out$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
                  conditional=FALSE, Ne=150000, ngen=100000, verbose=FALSE)
    d <- mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8))
    
    runid <- res
    output <- data.frame(id=runid, mu=d[1], nu=d[2], s=d[3], feature=fea)
    
    save(list=c("output", runid), file=paste0(outdir, "/", runid, ".RData"))
    write.table(output, paste0(outdir, "/", runid,".csv"), sep=",", row.names=FALSE, quote=FALSE)
    return(res)
}

####################################################################

### find overlap between genomic feature and COMET
### use findOverlaps from GenomicsRanges, return unqiue COMET blocks overlap with genomic features
get_overlap_sfs <- function(grf, grc, df){
    ## grf: GR feature object
    ## grc: GR Comet
    ## df: original data.frame for SFS
    
    tb <- findOverlaps(query=grf, subject=grc)
    tb <- as.matrix(tb)
    out <- as.data.frame(grc[unique(tb[,2])])
    df1 <- subset(df, cid %in% out$cid)
    sfs <- as.data.frame(table(df1$sfs))
    message(sprintf("[get_overlap_sfs]: [ %s ] sites of [ %s ] rows ", sum(sfs[,2]), nrow(sfs)))
    if(nrow(sfs) != 41){
        stop("[get_overlap_sfs]: not 41 sites!")
    }else{
        return(sfs)
    }
}
