
files <- list.files(path="profiling/2.Teo20_WGBS", pattern="0.E.", full.names=T)
file.rename(from=files, to=gsub("0.E.", "2.A", files))
