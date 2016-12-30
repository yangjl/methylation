### Jinliang Yang
### 10-12-2016
### COMET plot summary


library(ggplot2)
library(tidyr)

cg <- read.csv("cache/COMET_CG_sum_stat.csv")
cg$tot <- cg$high + cg$med + cg$low
cg$high1 <- cg$high/cg$tot
cg$med1 <- cg$med/cg$tot
cg$low1 <- cg$low/cg$tot


res2 <- cg %>% gather(key=level, value="ratio", 6:8)
res2$context <- toupper(res2$context)

res3 <- data.frame(level=c("high", "med", "low"), value=c(mean(cg$high1), mean(cg$med1), mean(cg$low1) ))
##### start to plot
p1 <- ggplot(res2, aes(x=level, y=ratio, fill=level)) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size=18, face="bold")) +
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0")) +
    scale_x_discrete(labels = c("High", "Med", "Low")) +
    ggtitle("COMET (BP)") + xlab("") + ylab("Ratio") + 
    guides(fill=FALSE)
#guides(colour=FALSE, linetype=FALSE)
p1
#########
pie(res3$value, labels=res3$level)

source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

pdf("graphs/stat_updated.pdf", width=16, height=5)
multiplot(p1, p2, cols=2)
dev.off()