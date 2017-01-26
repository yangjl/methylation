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
out$q <- gsub(".*_", "", out$id)

write.csv(out, "reports/rpkm_mean_pars.csv")


fit1 <- lm(nes ~ context + feature + q, data = out )
anova(fit1)


out <- read.csv("reports/rpkm_mean_pars.csv")

##########
files <- list.files(path="largedata/lcache", pattern="rpkm_var.*csv", full.names=TRUE)

out <- data.frame()
for(i in 1:length(files)){
    tb <-  read.csv(files[i])
    out <- rbind(out, tb)
}

out$nes <- 150000*out$s

out$context <- gsub("_.*", "", out$id)
out$q <- gsub(".*_", "", out$id)

write.csv(out, "reports/rpkm_var_pars.csv")


fit1 <- lm(nes ~ context + feature + q, data = out )
anova(fit1)


out <- read.csv("reports/rpkm_mean_pars.csv")

