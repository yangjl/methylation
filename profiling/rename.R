
files <- list.files(path="profiling/5.MCMC_scan/", pattern="4.H.", full.names=T)
file.rename(from=files, to=gsub("4.H.", "5.A", files))

files <- list.files(path="profiling/5.MCMC_scan/", pattern="5.A", full.names=T)
file.rename(from=files, to=gsub("5.A", "5.A.", files))
