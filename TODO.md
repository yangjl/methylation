## To Do List:
1. plot the genome-wide methylation differences.
2. identify DMR
3. correlate DMR with genomic features
4. Whether selection sweeps associated with differential methylation patterns: see new Tajima's D-like statistic

## Problems
1. fix the Comethylation boundary and missing data issue
2. IBD failure:
```
Module slurm/16.05.5 loaded
Module openmpi/2.0.1 loaded
Module JAVA 1.8 Loaded.
Exception in thread "main" java.lang.IllegalArgumentException:
8       174956473       8_174956473     C       G,*
7       1254    7_1254  T       G
        at ibd.IbdSegment.checkArguments(IbdSegment.java:85)
        at ibd.IbdSegment.<init>(IbdSegment.java:68)
        at main.WindowWriter.merge(WindowWriter.java:248)
        at main.WindowWriter.printIbd(WindowWriter.java:226)
        at main.WindowWriter.printIbd(WindowWriter.java:205)
        at main.Main.printOutput(Main.java:207)
        at main.Main.phaseData(Main.java:157)
        at main.Main.main(Main.java:115)
```


## QC report of WGBS

[1] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRA1_CTTGTA_R1_fastqc.html"  
[2] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRA1_CTTGTA_R2_fastqc.html"  
[3] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRA2_TGACCA_R1_fastqc.html"  
[4] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRA2_TGACCA_R2_fastqc.html"  
[5] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRA3_CCGTCC_R1_fastqc.html"  
[6] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRA3_CCGTCC_R2_fastqc.html"  
[7] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRB1_NoIndex_R1_fastqc.html" 
[8] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRB1_NoIndex_R2_fastqc.html" 
[9] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRB2_NoIndex_R1_fastqc.html"  
[10] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRB2_NoIndex_R2_fastqc.html"  
[11] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRB3_GTCCGC_R1_fastqc.html"   
[12] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRB3_GTCCGC_R2_fastqc.html"   
[13] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRC1_AGTTCC_R1_fastqc.html"  
[14] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRC1_AGTTCC_R2_fastqc.html"   
[15] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRC2_GCCAAT_R1_fastqc.html"  
[16] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRC2_GCCAAT_R2_fastqc.html"  
[17] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRC3_GTGAAA_R1_fastqc.html"  
[18] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRC3_GTGAAA_R2_fastqc.html"  
[19] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRD1_ATGTCA_R1_fastqc.html"  
[20] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRD1_ATGTCA_R2_fastqc.html"  
[21] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRD2_CAGATC_R1_fastqc.html"  
[22] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRD2_CAGATC_R2_fastqc.html"  
[23] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRD3_NoIndex_R1_fastqc.html" 
[24] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRD3_NoIndex_R2_fastqc.html" 
[25] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRE1_CCGTCC_R1_fastqc.html"  
[26] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRE1_CCGTCC_R2_fastqc.html"  
[27] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRE2_CTTGTA_R1_fastqc.html"  
[28] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRE2_CTTGTA_R2_fastqc.html"  
[29] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRF1_NoIndex_R1_fastqc.html"  
[30] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRF1_NoIndex_R2_fastqc.html" 
[31] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRF2_AGTCAA_R1_fastqc.html"  
[32] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRF2_AGTCAA_R2_fastqc.html"  
[33] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRG1_GTGAAA_R1_fastqc.html"  
[34] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRG1_GTGAAA_R2_fastqc.html"  
[35] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRG2_AGTTCC_R1_fastqc.html"  
[36] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRG2_AGTTCC_R2_fastqc.html"  
[37] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRH1_CGATGT_R1_fastqc.html" 
[38] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRH1_CGATGT_R2_fastqc.html"  
[39] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRH2_ATGTCA_R1_fastqc.html"  
[40] "http://htmlpreview.github.io/?https://github.com/yangjl/methylation/blob/master/reports/wgbs_qc_report/JRH2_ATGTCA_R2_fastqc.html" 
