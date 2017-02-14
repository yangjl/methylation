### Jinliang Yang
### 02-03-2017
### purpose: plot with the summarized methylation levels

##get command line args

cg <- read.csv("cache/CG_chr1_gff_feas.csv")


up1k <- subset(cg, feature %in% "up1k")
up1k$x <- gsub("V", "", up1k$var)
up1k$x <- as.numeric(as.character(up1k$x)) - 10


exon <- subset(cg, feature %in% "exon")
exon$x <- as.numeric(as.character(gsub("V", "", exon$var)))

intron <- subset(cg, feature %in% "intron")
intron$x <- gsub("V", "", intron$var)
intron$x <- as.numeric(as.character(intron$x)) + 10



down1k <- subset(cg, feature %in% "down1k")
down1k$x <- gsub("V", "", down1k$var)
down1k$x <- as.numeric(as.character(down1k$x)) + 20


ad <- rbind(up1k, exon, intron, down1k)


gene <- subset(cg, feature %in% "gene")
gene$x <- gsub("V", "", gene$var)
gene$x <- as.numeric(as.character(gene$x))

ad <- rbind(up1k, gene, down1k)
plot(ad$x, ad$mc)

ad$context <- "CG"

library(ggplot2)
library(reshape2)
source("lib/multiplot.R")

plot_eff <- function(outfile, getpdf){
    
    #######
    theme_set(theme_grey(base_size = 18)) 
    
    fsize=18
    p1 <- ggplot(ad, aes(x=x, y=mc, colour=factor(context)) )+
        labs(colour="context") +
        theme_bw() +
        xlab("GERP Score") +
        ylab("Additive Effect") +
        #scale_color_manual(values=cols) +
        #scale_linetype_manual(values=lty1) +
        guides(colour=FALSE, linetype=FALSE) +
        geom_smooth(span = 0.05) +
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize) )
    p1
    
}

########
p <- plot_eff(outfile="graphs/Fig2b_var.pdf", getpdf)
p
if(getpdf){
    pdf("graphs/Fig2b_var.pdf", width=8, height=4)
    p
    dev.off()
}