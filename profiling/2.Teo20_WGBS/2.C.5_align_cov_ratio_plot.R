### Jinliang Yang
### 09-06-2016

### CG, CHG and CHH
###180125000, 158277169, 624401016
library(ggplot2)
library(tidyr)

res1 <- read.csv("cache/wgbs_cov_09062016.csv")
res1$ratio <- with(res1, methylated/total)

res2 <- read.csv("cache/wgbs_ratio_09062016.csv")
res2 <- res2 %>% gather(key=context, value="ratio", 2:4)
res2$context <- toupper(res2$context)

##### start to plot
p1 <- ggplot(res1, aes(x=context, y=ratio, fill=context)) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size=18, face="bold")) +
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0")) +
    ggtitle("Coverage of cytosine sites") + xlab("") + ylab("covered / possible C sites") + 
    #ggtitle("Sequencing Depth") + xlab("") + ylab("Depth per cytosine site") + 
    guides(fill=FALSE)
#guides(colour=FALSE, linetype=FALSE)

p2 <- ggplot(res2, aes(x=context, y=ratio, fill=context)) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size=18, face="bold")) +
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0")) +
    ggtitle("Sequencing Depth") + xlab("") + ylab("Depth per cytosine site") + 
    #ggtitle("Methylation Ratio") + xlab("") + ylab("Mean C/CT Ratio") + 
    guides(fill=FALSE)

p3 <- ggplot(subset(resl, type=="tot"), aes(x=context, y=value, fill=context)) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.title = element_text(size=18, face="bold")) +
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0")) +
    ggtitle("Coverage of cytosine sites") + xlab("") + ylab("covered / possible C sites") + 
    guides(fill=FALSE)
#guides(colour=FALSE, linetype=FALSE)

source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

pdf("graphs/stat_updated.pdf", width=16, height=5)
multiplot(p1, p2, cols=2)
dev.off()
