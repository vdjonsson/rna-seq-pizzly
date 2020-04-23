log <- file(snakemake@output[[1]], open="wt")
sink(log)
sink(log, type="message")

library(rhdf5)
library(tidyverse)

samples <- read_tsv(snakemake@input[["samples"]], na = "", col_names = TRUE) %>%
            # make everything except the sample name and path string a factor
            mutate_at(  vars(-sample, -path),
                        list(~factor(.))
                        )


s2c <- samples[c("sample", "condition", "path")]
apply(s2c,1,function(f){ dh5 <- try(dim(h5ls(paste0(f["path"],"/abundance.h5")))[1]); if(dh5!=115){ dh5<-"ERROR" }; return(paste0(f["path"]," -> ",dh5)) })
