### Jinliang Yang
### Jan 10th, 2017


files <- list.files(path="largedata/lcache", pattern="rpkm_mean.*csv", full.names=TRUE)

out <- data.frame()
for(i in 1:length(files)){
   tb <-  read.csv(files[i])
   out <- rbind(out, tb)
}

out$nes <- 150000*out$s

out$context <- gsub("_.*", "", out$id)
out$q <- paste("q", gsub(".*_", "", out$id), sep="")

write.csv(out, "reports/gbody_popgen_pars.csv")

out <- read.csv("reports/popgen_pars.csv")

##########
files <- list.files(path="largedata/lcache", pattern="rpkm_var.*csv", full.names=TRUE)

out <- data.frame()
for(i in 1:length(files)){
    tb <-  read.csv(files[i])
    out <- rbind(out, tb)
}

out$nes <- 150000*out$s

out$context <- gsub("_.*", "", out$id)
out$q <- paste("q", gsub(".*_", "", out$id), sep="")

write.csv(out, "reports/gbody_popgen_pars.csv")

out <- read.csv("reports/popgen_pars.csv")
