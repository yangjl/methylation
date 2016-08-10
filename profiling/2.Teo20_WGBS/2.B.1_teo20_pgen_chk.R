### Jinliang Yang
### August 8th, 2016


### filter to biallelic loci
# bcftools view JRI20_filtered_snps_annot.bcf.gz -m2 -M2 -v snps -Oz -o JRI20_bi_snps_annot.vcf.gz
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%IUPACGT]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.txt
# bcftools query -f 'chr\tpos\tref\talt[\t%SAMPLE]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.header

library("data.table")
h <- read.table("largedata/gatk_vcf/JRI20_bi_snps_annot.header", header=T)
snp <- fread("largedata/gatk_vcf/JRI20_bi_snps_annot.txt", nrows=10000)
names(snp) <- names(h)

library("Biostrings")
fa <- readDNAStringSet(filepath = "~/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa", format="fasta")
bck <- fa
uniqueLetters(fa)
#[1] "A" "C" "G" "T" "S" "Y" "N"
alphabetFrequency(fa)
# length of each chromosome
width(fa)
names(fa)


sub_rules <- data.frame(from=c("M", "Y", "R", "K"), to=c("C", "C", "G", "T"))



sub_bp <- function(dt, ){
    
    
    
    
    
    return(dt)
}







### replace bp

chr1

res <- replaceAt(x=fa[13], at=1:10, value="M")
subseq(fa[13], start=1, end=10)











test <- sort(fa)
### edit genome to remove 
test <- fa[-1:-3]




letterFrequency(fa, letters="ATCG")
names(fa)

writeXStringSet()


uniqueLetters(fa)
# "*" "A" "C" "G" "K" "M" "N" "R" "S" "T" "W" "Y"

fa[12] %in% "MM"

subseq, subseq<- extractAt, replaceAt


vmatchPattern(pattern="*", fa[4])



get_bp <- function(pos0=987, chr=12){
    print(subseq(fa[chr], start=pos0, end=pos0))
}
get_bp(pos0=92673)
