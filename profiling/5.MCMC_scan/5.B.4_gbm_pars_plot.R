







outall <- read.csv("reports/gbody_popgen_pars.csv")
fit1 <- lm(nes ~ context + feature + q, data = outall )
anova(fit1)

fit2 <- lm(mu ~ context + feature + q, data = outall )
anova(fit2)

fit3 <- lm(nu ~ context + feature + q, data = outall )
anova(fit3)

fit4 <- lm(nu ~ mu + context + feature + q, data = outall ) 
anova(fit4)

library(ggplot2)
library(cowplot)

plot_fig3ab <- function(outfile, getpdf){
    outall <- read.csv("reports/gbody_popgen_pars.csv")
    
    out <- outall
    fsize=16
    p1 <- ggplot(out, aes(x=factor(feature), y= log10(mu),
                          fill=factor(context, levels=c("CG", "CHG"), labels=c("CG", "CHG")))) + 
        geom_boxplot() +
        facet_wrap(~ q, ncol = 1) + 
        xlab("") +
        guides(fill=FALSE) +
        ylim(c(-7, -5)) +
        ylab("log10(mu)") +
        ggtitle("") + theme_bw() +
        labs(fill="Context") +
        theme(axis.text = element_text( angle = 0, size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    p2 <- ggplot(out, aes(x=factor(feature), y= log10(nu),
                          fill=factor(context, levels=c("CG", "CHG"), labels=c("CG", "CHG")))) + 
        geom_boxplot() +
        facet_wrap(~ q, ncol = 1) + 
        guides(fill=FALSE) +
        xlab("") +
        ylim(c(-7, -5)) +
        ylab("log10(nu)") +
        ggtitle("") + theme_bw() +
        labs(fill="Context") +
        theme(axis.text = element_text(angle = 0, size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    p3 <- ggplot(out, aes(x=factor(feature), y= nes,
                          fill=factor(context, levels=c("CG", "CHG"), labels=c("CG", "CHG")))) + 
        geom_boxplot() +
        facet_wrap(~ q, ncol = 1) + 
        xlab("") +
        #ylim(c(0,1)) +
        ylab("Ne*s") +
        ggtitle("") + theme_bw() +
        labs(fill="Context") +
        theme(axis.text = element_text(angle = 0, size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    #PLOT
    plot_grid(p1, p2, p3, ncol=3, rel_widths=c(1,1,1.3), align="h", labels =c("A", "B", "C"))
}
######
plot_fig3ab(outfile="graphs/Fig_post_var.pdf", getpdf)
