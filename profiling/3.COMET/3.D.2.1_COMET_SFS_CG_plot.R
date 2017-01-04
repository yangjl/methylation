### Jinliang Yang
### 12-07-2016
### plot sfs with multiple features

#####>>> read from cache/
library("data.table")
library("tidyr")
library("ggplot2")
library(cowplot)

cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###################################################################################################
### SFS and bp
df <- read.csv("cache/SFS_comet_blocks_CG.csv")
dt <- as.data.table(df)
tab1 <- dt[, .(bp = sum(length)), by= sfs] 
tab2 <- data.frame(table(df$sfs))

tab1 <- as.data.frame(tab1)
plot(tab1$sfs, tab1$bp, type="h")
tab <- merge(tab1, tab2, by.x="sfs", by.y="Var1")

tab$bp <- tab$bp/sum(tab$bp)
tab$Freq <- tab$Freq/sum(tab$Freq)

out <- gather(tab, "type", "value", 2:3)

s1 <- ggplot(out, aes(sfs, value, fill=type)) + geom_bar(stat="identity", position = "dodge") +
    #theme_bw() +
    xlab("") +
    ylab("Frequency") +
    #scale_fill_manual(values=c("#ff0000", "#008080", "#003366"),
    #                  name="", labels=c("Class I", "Class II and III", "Other TE")) +
    theme(legend.position="top", axis.text=element_text(size=15), axis.title=element_text(size=15) ) +
    scale_fill_discrete(name="", labels=c("Base-pairs", "COMET frequency"))
#########
s1

###################################################################################################
### SFS and genic features
gen <- read.csv("cache/SFS_comet_features.csv")
gen$exon <- gen$exon/sum(gen$exon)
gen$intron <- gen$intron/sum(gen$intron)
gen$onek <- gen$onek/sum(gen$onek)

out2 <- gather(gen[, c("site", "onek", "exon", "intron")], "type", "value", 2:4)
s2 <- ggplot(out2, aes(site, value, fill=type)) +
    geom_bar(stat="identity", position = "dodge") +
    #theme_bw() +
    xlab("") +
    ylab("Frequency") +
    #scale_fill_manual(values=c("#ff0000", "#008080", "#003366"),
    #                  name="", labels=c("Exon", "Intron", "1k upstream")) +
    scale_fill_discrete(name="", labels=c("Exon", "Intron", "1k upstream")) +
    theme(legend.position="top", axis.text=element_text(size=15), axis.title=element_text(size=15) )
#########
s2

###################################################################################################
### SFS and transposons
tab4 <- read.csv("cache/SFS_comet_features.csv")
tab4$c1 <- tab4$c1/sum(tab4$c1)
tab4$c2 <- tab4$c2/sum(tab4$c2)
tab4$c3 <- tab4$c3/sum(tab4$c3)

out3 <- gather(tab4[, c("site", "c1", "c2", "c3")], "type", "value", 2:4)
s3 <- ggplot(out3, aes(site, value, fill=type)) +
    geom_bar(stat="identity", position = "dodge") +
    #theme_bw() +
    xlab("") +
    ylab("Frequency") +
    #scale_fill_manual(values=c("#ff0000", "#008080", "#003366"),
    #                  name="", labels=c("Class I", "Class II and III", "Other TE")) +
    scale_fill_discrete(name="", labels=c("Class I", "Class II and III", "Other TE")) +
    theme(legend.position="top", axis.text=element_text(size=15), axis.title=element_text(size=15) ) 
#########
s3


s4 <- ggplot(df, aes(factor(sfs), log(length))) + geom_boxplot() +
    #theme_bw() +
    xlab("") +
    ylab("COMET log(bp)") +
    #scale_fill_manual(values=c("#ff0000", "#008080", "#003366"),
    #                  name="", labels=c("Class I", "Class II and III", "Other TE")) +
    scale_x_discrete(breaks=c(0, 10, 20, 30, 40)) +
    theme(legend.position="top", axis.text=element_text(size=15), axis.title=element_text(size=15) )

#########
s4

###################################################################################################
### SFS and transposons
### cowplot: combined
p <- plot_grid(s1, s2, s3, s4, ncol=1, rel_heights=c(1,1,1,1), align="v")

pdf("graphs/Figure_sfs_cg.pdf", width=6, height=12)
ggdraw(p) + draw_plot_label(c("A", "B", "C", "D" ), c(0, 0, 0, 0), c(1, 3/4, 2/4, 1/4), size=15)
dev.off()





plot_grid(s1,ncol=1, align="v")


