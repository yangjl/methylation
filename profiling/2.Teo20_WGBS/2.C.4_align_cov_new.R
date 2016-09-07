### Jinliang Yang
### plot the stat of the teo20 methylation data


### write a python script and run through all the samples
library("farmeR")

files <- list.files(path="largedata/wgbs_align", pattern="pe.CX_report.txt", full.names=TRUE)

df <- data.frame(input=files, output=files)
df$output <- gsub("align", "smoothed", df$output)
df$output <- gsub(".CX_report.txt", "", df$output)

for(i in 1:19){
    shid <- paste0("slurm-script/split_", i, ".sh")
    command <- paste("python lib/splitcontext.py -i", df$input[i], "-o", df$output[i])
    cat(command, file=shid, sep="\n", append=FALSE)
}
shcode <- paste("sh slurm-script/split_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

set_array_job(shid="slurm-script/run_split.sh", shcode=shcode,
              arrayjobs="1-19", wd=NULL, jobid="split",
              email= "yangjl0930@gmail.com",
              run = c(TRUE, "bigmemm", 1))


### collect results and plot it

get_cov_ratio <- function(){
    ### get coverage
    files <- list.files(path="largedata/wgbs_smoothed", pattern="stat", full.names=TRUE)
    out1 <- data.frame()
    for(i in 1:20){
        f <- read.table(files[i], header=TRUE)
        f$file <- files[i]
        out1 <- rbind(out1, f)
    }
    write.table(out1, "cache/wgbs_cov_09062016.csv", sep=",", row.names=FALSE, quote=FALSE)
    
    
    ### get ratio of CG, CHG and CHH
    fcg <- list.files(path="largedata/wgbs_smoothed", pattern="cg", full.names=TRUE)
    fchg <- list.files(path="largedata/wgbs_smoothed", pattern="chg", full.names=TRUE)
    fchh <- list.files(path="largedata/wgbs_smoothed", pattern="chh", full.names=TRUE)
    
    out2 <- data.frame()
    for(i in 1:20){
        f1 <- fread(fcg[i])
        r1 <- (f1[,sum(V4)] + f1[,sum(V5)])/nrow(f1)
        
        f2 <- fread(fchg[i])
        r2 <- (f2[,sum(V4)] + f2[,sum(V5)])/nrow(f2)
        
        f3 <- fread(fchh[i])
        r3 <- (f3[,sum(V4)] + f3[,sum(V5)])/nrow(f3)
        tem <- data.frame(file=fcg[i], cg=r1, chg=r2, chh=r3)
        out2 <- rbind(out2, tem)
    }
    write.table(out2, "cache/wgbs_ratio_09062016.csv", sep=",", row.names=FALSE, quote=FALSE)
}

#######
get_cov_ratio()




