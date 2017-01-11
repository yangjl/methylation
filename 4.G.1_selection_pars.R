### Jinliang Yang
### Jan 10th, 2017


files <- list.files(path="largedata/lcache", pattern="^C.*", full.names=TRUE)

out <- data.frame()
for(i in 1:length(files)){
   tb <-  read.csv(files[i])
   out <- rbind(out, tb)
}

out$nes <- 150000*out$s

out$context <- gsub("_.*", "", out$id)
out$q <- paste("q", gsub(".*_", "", out$id), sep="")


fit1 <- lm(mu ~ context + feature + q, data = out )
anova(fit1)

fit2 <- lm(nu ~ context + feature + q, data = out )
anova(fit2)

fit3 <- lm(s ~ context + feature + q, data = out )

anova(fit3)


fit4 <- lm(nu ~ mu + context + feature + q, data = out )
anova(fit4)

fit5 <- lm(mu ~ nu + context + feature + q, data = out )
anova(fit5)

fit6 <- lm(s ~ mu + nu + context + q + feature, data = out )
anova(fit6)

