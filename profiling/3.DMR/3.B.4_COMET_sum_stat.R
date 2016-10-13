### Jinliang Yang
### 10-12-2016
### COMET summary statistics

cat_comet <- function(){
    dirs <- dir(path="largedata/COMET", pattern="JR", full.names=TRUE, recursive=FALSE)
    
    for(i in 1:length(dirs)){
        
        files <- list.files(path=paste0(dirs[i], "/COMETs"), pattern="txt", full.names=TRUE)
        
        out <- data.frame()
        sid <- gsub(".*/", "", dirs[i])
        for(j in 1:length(files)){
            chr <- read.table(files[j], header=TRUE)
            out <- rbind(out, chr)
        }
        out$level <- "med"
        out[out$meth >= 0.66, ]$level <- "high"
        out[out$meth <= 0.33, ]$level <- "low"
        
        message(sprintf("###>>> writing [ %s ] ...", sid))
        write.table(out, paste0("largedata/COMET/CG_COMET/", sid, "_CG_COMET.csv"), 
                    sep=",", row.names=FALSE, quote=FALSE)
    }
}

##########
cat_comet()






