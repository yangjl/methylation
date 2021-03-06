### Jinliang Yang
### Convert the BSMAP to BED12 format

indir= "/group/jrigrp4/BS_teo20/BSMAP_round2"
outdir= "/home/jolyang/Documents/Github/methylation/largedata/vcf_files/"
files=list.files(path=indir, pattern="txt$", full.names = TRUE)
inputdf <- data.frame(bsmap=files, out=paste0(outdir, gsub(".*/|txt", "", files), "vcf") )

library(farmeR)
run_bs2vcf(inputdf, email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 1))
## bs2vcf -p /group/jrigrp4/BS_teo20/BSMAP_output -i test_methratio.txt -o test_methratio.vcf

run_bgzip <- function(indir="/group/jrigrp4/BS_teo20/BSMAP_output"){
    
    files=list.files(path=indir, pattern="methratio.vcf$")
    for(i in 1:length(files)){
        shid <- paste0("slurm-script/run_bgzip_", i, ".sh")
        #out <- gsub(".*/", "", out)
        #outfile <- paste0(outdir, "/", out)
        cmd1 <- paste0("cd ", indir)
        cmd2 <- paste0("bgzip ", files[i],  "; tabix -p vcf ", files[i], ".gz")
        cat(c(cmd1, cmd2), file=shid, sep="\n", append=FALSE)
    }
    set_array_job(shid = "slurm-script/run_bgzip.sh",
                  shcode = "sh slurm-script/run_bgzip_$SLURM_ARRAY_TASK_ID.sh", 
                  arrayjobs = "1-20", wd = NULL,
                  jobid = "bgzip", email = "yangjl0930@gmail.com", 
                  runinfo = c(FALSE, "bigmemm", "1"))
}

#bgzip test_methratio.vcf; tabix -p vcf test_methratio.vcf.gz
run_bgzip(indir="~/Documents/Github/methylation/largedata/vcf_files")
###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8196 --ntasks=1 slurm-script/run_bgzip.sh

system("ls | grep 'vcf.gz' > vcflist.txt")
#bcftools merge -l vcflist.txt -o teo12_methratio.vcf -O b

library(farmeR)
shcode <- c("cd ~/Documents/Github/methylation/largedata/vcf_files", 
            "bcftools merge -l vcflist.txt -o teo20_methratio.bcf -O b --force-samples",
            "bcftools index teo20_methratio.bcf")
set_farm_job(slurmsh="slurm-script/run_bcf_merge.sh", shcode=shcode,
           wd=NULL, jobid="bcfmerge", email="yangjl0930@gmail.com",
           runinfo=c(TRUE, "bigmemh", "2"))
###>>> sbatch -p bigmemh --mem 16392 --ntasks=2 slurm-script/run_bcf_merge.sh

###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks=2 mem 16000 slurm-script/run_bcf_merge.sh

# bcftools index bcftools index teo12_methratio.vcf.gz


#bcftools query -f '%ID\t%CO[\t%RA]\n' -r 1:0-10000 teo20_methratio.bcf -o test_RA.txt
library(farmeR)
shcode <- c("cd ~/Documents/Github/methylation/largedata/vcf_files",
            "bcftools query -f \'%ID\\t%CO[\\t%RA]\\n\' teo20_methratio.bcf -o teo20_RA.txt")
set_farm_job(slurmsh="slurm-script/run_bcfra.sh", shcode=shcode,
             wd=NULL, jobid="bcfra", email="yangjl0930@gmail.com",
             runinfo=c(FALSE, "bigmemh", "2"))