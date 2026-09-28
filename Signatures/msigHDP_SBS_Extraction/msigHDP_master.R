#!/usr/bin/env Rscript

library(mSigHdp)
library(ICAMS)
library(cosmicsig)

#note: the input file require the trinucleotides as format ACTA, meaning A[C>A]T 

input_data <- read.table("./01_Input/msigHDP_input.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)

input_data <- ICAMS::as.catalog(input_data, catalog.type = "counts")

results <- mSigHdp::RunHdpxParallel(
  input.catalog        = round(input_data),
  out.dir              = "./02_Output",
  num.child.process    = 20, 
  CPU.cores            = 20,
  seedNumber           = 42,
  K.guess              = 3, #3 found in the original study, increased for large cohort size
  burnin               = 5000,
  burnin.multiplier    = 20,
  post.n               = 200, 
  post.space           = 100, 
  multi.types          = FALSE,
  overwrite            = TRUE,
  gamma.alpha          = 1,
  gamma.beta           = 20, 
  high.confidence.prop = 0.7,
  checkpoint           = TRUE,
  verbose              = TRUE) 

save(results, file = "./02_output/multi.chains.Rdata")
