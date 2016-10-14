### Jinliang Yang
### 10-12-2016
### Chop COMET into shared blocks


comet2blocks <- function(files, chri=10, cutoff=c(0.33, 0.66)){
    
    cometls <- list()
    for(i in 1:length(files)){
        df <- fread(files[i], data.table=FALSE)
        df <- subset(df, chr == chri)
        
        message(sprintf("[comet2blocks] recalculating low, med, high, using cutoff [ <=%s and >=%s] ... ",
                        cutoff[1], cutoff[2]))
        df$level <- "med"
        df[df$meth <= cutoff[1], ]$level <- "low"
        df[df$meth >= cutoff[2], ]$level <- "high"
        
        sid <- gsub(".*/|_.*", "", files[i])
        message(sprintf("[comet2blocks] condense [sample %s] for [chr %s] ... ",
                        sid, chri))
        cometls[[sid]] <- condense(df)
    }
    
    #message(sprintf("[comet2blocks] chop into shared blocks ... ")
    res <- chop2blocks(cometls, chri)
    res <- apply(res, 2, as.character)
    res[res=="high"] <- 2
    res[res=="med"] <- 1
    res[res=="low"] <- 0
    res[res=="a"] <- -9
    return(res)
}



## condense COMET by chr
condense <- function(df){
    df <- df[order(df$start), ]

    df$level <- as.character(df$level)
    df$level0 <- c("no", df$level[-nrow(df)])
    df$eval <- 1
    df[df$level == df$level0, ]$eval <- 0
    df$block <- cumsum(df$eval)
    
    out <- ddply(df, .(block, level), summarise,
                 start=min(start),
                 end=max(stop),
                 meth=mean(meth))
    return(out)
}

## chop into shared blocks
chop2blocks <- function(cometls, chri){
    
    ##### determine the break points
    out <- lapply(1:length(cometls), function(x){
        return(unique(c(cometls[[x]]$start, cometls[[x]]$end)))
    })
    bp <- c()
    for(i in 1:length(out)){
        bp <- c(bp, out[[i]])
    }
    bp <- sort(unique(bp))
    message(sprintf("[chop2blocks] chop [chr%i] into [ %s ] shared blocks ... ", chri, length(bp)))
    
    ### use genomicranges to find overlaps
    gr <- GRanges(seqnames=Rle(paste0("chr", chri)),
                  ranges=IRanges(start=bp[-length(bp)], end=bp[-1]-1),
                  strand = Rle(strand("+")) )
    
    out <- data.frame(chr=seqnames(mygr), bid=paste(start(mygr), end(mygr), sep="_"))
    for(j in 1:length(cometls)){
        message(sprintf("[chop2blocks] working on [ sample: %s] ... ", names(cometls[j])))
        ### use genomicranges to find overlaps
        gr1 <- GRanges(seqnames=paste0("chr", chri),
                       ranges=IRanges(start=cometls[[j]]$start, end=cometls[[j]]$end-1),
                       strand = Rle(strand("+")),
                       score = cometls[[j]]$meth,
                       level = cometls[[j]]$level)
        mygr <- gr
        idx <- findOverlaps(query=gr1, subject=mygr, select="all")
        #mygr$score <- -9
        #mygr[subjectHits(idx), ]$score <- gr1[queryHits(idx), ]$score
        mygr$level <- "a"
        mygr[subjectHits(idx), ]$level <- gr1[queryHits(idx), ]$level
        
        mx <- data.frame( bid=paste(start(mygr), end(mygr), sep="_"), level=mcols(mygr)$level)
        names(mx)[2] <- names(cometls[j])
        out <- merge(out, mx, by="bid", sort=FALSE)
    }
    return(out)
}



