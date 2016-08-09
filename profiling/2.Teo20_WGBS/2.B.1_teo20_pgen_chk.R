### Jinliang Yang
### August 8th, 2016

library("Biostrings")
fa <- readBStringSet(filepath = "largedata/wgbs_pgen/JRIAL2A/JRIAL2A.fa", format="fasta")
alphabetFrequency(fa, letters="ACGT")
letterFrequency(fa, letters="ATCG")
names(fa)

uniqueLetters(fa)
# "*" "A" "C" "G" "K" "M" "N" "R" "S" "T" "W" "Y"

fa[12] %in% "MM"

subseq, subseq<- extractAt, replaceAt


vmatchPattern(pattern="*", fa[4])



get_bp <- function(pos0=987, chr=12){
    print(subseq(fa[chr], start=pos0, end=pos0))
}
get_bp(pos0=92673)
