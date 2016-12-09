### Jinliang Yang
### 10-14-2016
### run MCMC for comet on chr10

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

source("lib/mplots.R")
source("lib/mcmcbc.R")
library(data.table)


run_mcmc <- function(myi, b){
    output <- data.frame()
    if(myi > 10) stop("[run_mcmc] error! [myi=%s] should not > 10!", myi)
    
    message(sprintf("[run_mcmc] starting for [ %s ] analysis ...", myi))
    ################################################################
    # SFS computed from COMET with varing length quantiles
    if(myi %in% 1:5){
        rdat_id <- paste0("pres", myi)
        if(myi == 1){
            ### using all data
            b0 <- b
            temp <- data.frame(type="CG", region="GW", size="all", bin="COMET", rdat=rdat_id)
        }
        if(myi == 2){
            b0 <- subset(b, length <= quantile(b$length)[2])
            temp <- data.frame(type="CG", region="GW", size="0-0.25", bin="COMET", rdat=rdat_id)
        }
        
        if(myi == 3){
            b0 <- subset(b, length < quantile(b$length)[2] & length <= quantile(b$length)[3])
            temp <- data.frame(type="CG", region="GW", size="0.25-0.5", bin="COMET", rdat=rdat_id)
        }
        
        if(myi == 4){
            b0 <- subset(b, length < quantile(b$length)[3] & length <= quantile(b$length)[4])
            temp <- data.frame(type="CG", region="GW", size="0.5-0.75", bin="COMET", rdat=rdat_id)
        }
        
        if(myi == 5){
            b0 <- subset(b, length > quantile(b$length)[4])
            temp <- data.frame(type="CG", region="GW", size="0.75-1", bin="COMET", rdat=rdat_id)
        }
        
        out <- table(b0$sfs)
        sfs <- data.frame(out)
        # If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
        res <- MCMCBC(my_sfs=sfs$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
                      conditional=FALSE, Ne=150000, ngen=100000, verbose=FALSE)
        d <- mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8))
        
        rdat_id <- res
        temp$mu <- d[1]
        temp$nu <- d[2]
        temp$s <- d[3]
        output <- rbind(output, temp)
        save(list=c("output", output$rdat), file=paste0("largedata/lcache/", rdat_id, "_mcmc_CG_GW.RData"))
        write.table(output, paste0("largedata/lcache/", rdat_id,"_CG_GW.csv"), sep=",", row.names=FALSE, quote=FALSE)
    }
    
    ################################################################
    # SFS caculated from bp
    if(myi %in% 6:10){
        dt <- as.data.table(b)
        rdat_id <- paste0("pres", myi)
        if(myi == 6){
            ### using all data
            b0 <- dt
            temp <- data.frame(type="CG", region="GW", size="all", bin="bp", rdat=rdat_id)
        }
        if(myi == 7){
            b0 <- subset(dt, length <= quantile(b$length)[2])
            temp <- data.frame(type="CG", region="GW", size="0-0.25", bin="bp", rdat=rdat_id)
        }
        
        if(myi == 8){
            b0 <- subset(dt, length < quantile(b$length)[2] & length <= quantile(b$length)[3])
            temp <- data.frame(type="CG", region="GW", size="0.25-0.5", bin="bp", rdat=rdat_id)
        }
        
        if(myi == 9){
            b0 <- subset(dt, length < quantile(b$length)[3] & length <= quantile(b$length)[4])
            temp <- data.frame(type="CG", region="GW", size="0.5-0.75", bin="bp", rdat=rdat_id)
        }
        if(myi == 10){
            b0 <- subset(b, length > quantile(b$length)[4])
            temp <- data.frame(type="CG", region="GW", size="0.75-1", bin="bp", rdat=rdat_id)
        }
        tab <- b0[, .(bp = sum(length)), by= sfs]
        tab <- tab[order(tab$sfs),]
        # If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
        res <- MCMCBC(my_sfs=tab$bp, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
                      conditional=FALSE, Ne=150000, ngen=100000, verbose=FALSE)
        d <- mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8))
        
        rdat_id <- res
        temp$mu <- d[1]
        temp$nu <- d[2]
        temp$s <- d[3]
        output <- rbind(output, temp)
        save(list=c("output", output$rdat), file=paste0("largedata/lcache/", rdat_id, "_mcmc_CG_GW.RData"))
        write.table(output, paste0("largedata/lcache/", rdat_id,"_CG_GW.csv"), sep=",", row.names=FALSE, quote=FALSE)
    }
    
    
    message("[run_mcmc] done!")
    
}
    
###### main codes:    
b <- read.csv("cache/SFS_comet_blocks_CG.csv")
#hist(log10(b$length), xlab="COMET Length log10(bp)", main="CG COMETs", col="#cdc0b0")
#abline(v=quantile(log10(b$length))[2:4], lwd=3, lty=2 )
# 0%   25%   50%   75%  100% 
# 1    96   309   784 76180 

run_mcmc(myi=JOBID, b)






