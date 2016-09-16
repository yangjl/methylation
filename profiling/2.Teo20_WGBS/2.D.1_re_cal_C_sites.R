### Jinliang Yang
### 09-13-2016
### purpose: re-calibrate the C/T and G/A heterozygote sites

library(data.table)

files <- list.files(path="largedata/wgbs_smoothed", pattern="cg$", full.names = TRUE)


re_calC <- function(){
    ### First, SNP table:
    h <- read.table("largedata/gatk_vcf/JRI20_bi_snps_annot.header", header=T)
    ### I do not need to merge because the orders are exactly the same.
    idtab <- read.csv("data/teo20_ids.csv")
    n <- c("chr", "pos", "ref", "alt", paste0("JR", idtab$plate))
    snpdt <- fread("largedata/gatk_vcf/JRI20_bi_snps_annot.txt")
    names(snpdt) <- n
    
    ### loop through genotype and then chr
    for(i in 1:length(files)){
        dt <- fread(files[i])
        
        for(j in 1:10){
            sid <- gsub(".*\\/|_pe.*", "", files[i])
            message(sprintf("###>>> re-cal sample [ %s, ID=%s ] chr [ %s ]", i, sid, j))
            chr <- dt[V1 == j]
            chr[, snpid := paste(V1, V2, sep="_")]
            
            
            sub <- snpdt[, c("chr", "pos", "ref", "alt", sid),  with=FALSE]
            sub <- sub[chr == 10]
            names(sub)[5] <- "SAMPLE"
            sub <- sub[SAMPLE == "Y" | SAMPLE == "R"]
            #sub$snpid <- paste(sub$chr, sub$pos, sep="_")
            
            sub[, snpid := paste(chr, pos, sep="_")]
            
            out <- merge(chr, sub, by.x="snpid", by.y="snpid", all.x=TRUE)
            out[, tot := V4 + V5]
            out[!is.na(SAMPLE), tot := ceiling((V4 + V5)/2)]
            out[V4 > tot & SAMPLE == "Y", tot := V4]
            
        }
        
    }
}





