
library(ggplot2)
source("lib/multiplot.R")

plot_fig3ab <- function(outfile, getpdf){
    out <- read.csv("reports/gbody_popgen_pars.csv")
    
    fsize=16
    p1 <- ggplot(out, aes(x=factor(feature), y= -log10(mu),
                           fill=factor(context, levels=c("CG", "CHG"), labels=c("CG", "CHG")))) + 
        geom_bar(position=position_dodge(), stat="identity") +
        xlab("") +
        #ylim(c(0,1)) +
        ylab("mu") +
        ggtitle("") + theme_bw() +
        labs(fill="Context") +
        theme(axis.text = element_text(angle = 90, hjust = 1, size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    p2 <- ggplot(out, aes(x=factor(feature), y= -log10(nu),
                          fill=factor(context, levels=c("CG", "CHG"), labels=c("CG", "CHG")))) + 
        geom_bar(position=position_dodge(), stat="identity") +
        xlab("") +
        #ylim(c(0,1)) +
        ylab("nu") +
        ggtitle("") + theme_bw() +
        labs(fill="Context") +
        theme(axis.text = element_text(angle = 90, hjust = 1, size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    p3 <- ggplot(out, aes(x=factor(feature), y= nes,
                          fill=factor(context, levels=c("CG", "CHG"), labels=c("CG", "CHG")))) + 
        geom_bar(position=position_dodge(), stat="identity") +
        xlab("") +
        #ylim(c(0,1)) +
        ylab("Ne*s") +
        ggtitle("") + theme_bw() +
        labs(fill="Context") +
        theme(axis.text = element_text(angle = 90, hjust = 1, size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    
}
######
plot_fig3ab(outfile="graphs/Fig_post_var.pdf", getpdf)
