### Jinliang Yang
### 10-16-2016
### compute IBD among 20 teo lines

library(farmeR)

shcode = c("module load java", "cd largedata/gatk_vcf/",
           paste("java -Xmx80g -jar ~/bin/beagle.22Feb16.8ef.jar",
                 "gt=JRI20_filtered_snps_annot.vcf.gz ibd=true out=JRI20_ibd"))

set_array_job(shid = "slurm-script/run_beagle.sh",
              shcode = shcode, arrayjobs = "1", wd = NULL,
              jobid = "beagle", email = "yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 10, "80G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 80G --ntasks=10 --time=48:00:00 slurm-script/run_beagle.sh