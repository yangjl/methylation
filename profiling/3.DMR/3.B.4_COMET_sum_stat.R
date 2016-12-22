### Jinliang Yang
### 10-12-2016
### COMET summary statistics

cat_comet <- function(path="largedata/COMET", outdir="largedata/COMET/CG_COMET/"){
    dirs <- dir(path=path, pattern="JR", full.names=TRUE, recursive=FALSE)
    ### checking the status
    for(i in 1:length(dirs)){
        files <- list.files(path=paste0(dirs[i], "/COMETs"), pattern="txt", full.names=TRUE)
        if(length(files) < 10){
            message(sprintf("###>>>> [ %s ] not finished smoothing!", dirs[i]))
        }
    }
    
    for(i in 1:length(dirs)){
        
        files <- list.files(path=paste0(dirs[i], "/COMETs"), pattern="txt", full.names=TRUE)
        
        out <- data.frame()
        sid <- gsub(".*/", "", dirs[i])
        for(j in 1:length(files)){
            chr <- fread(files[j], header=TRUE, data.table=FALSE)
            out <- rbind(out, chr)
        }
        out$level <- "med"
        out[out$meth >= 0.66, ]$level <- "high"
        out[out$meth <= 0.33, ]$level <- "low"
        
        message(sprintf("###>>> writing [ %s ], [ %s ] ...", sid, nrow(out)))
        write.table(out, paste0(outdir, sid, "_COMET.csv"), 
                    sep=",", row.names=FALSE, quote=FALSE)
    }
}

library("data.table")
##########
### CG
cat_comet(path="largedata/COMET", outdir="largedata/COMET/CG_COMET/")
### CHG
cat_comet(path="largedata/COMET_CHG", outdir="largedata/COMET_CHG/CHG_COMET/")




############################################################################
get_sum <- function(path="largedata/COMET/CG_COMET"){
    files <- list.files(path=path, pattern="csv", full.names=T)
    out <- data.frame()
    for(i in 1:length(files)){
        df <- read.csv(files[i])
        sid <- gsub(".*/|_.*", "", files[i])
        message(sprintf("###>>> Processing sample: [ %s ] ...", sid))
        tem <- data.frame(id=sid, high=sum(subset(df, level %in% "high")$size),
                          med=sum(subset(df, level %in% "med")$size),
                          low=sum(subset(df, level %in% "low")$size) )
        
        out <- rbind(out, tem)
    }
    return(out)
}

##########
## CG
res <- get_sum(path="largedata/COMET/CG_COMET")
write.table(res, "cache/COMET_CG_sum_stat.csv", sep=",", row.names=FALSE, quote=FALSE)

## CHG
res <- get_sum(path="largedata/COMET_CHG/CHG_COMET/")
write.table(res, "cache/COMET_CHG_sum_stat.csv", sep=",", row.names=FALSE, quote=FALSE)


stat <- read.csv("cache/COMET_CG_sum_stat.csv")
