## This script splits the large data into 5 equal parts
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(Matrix)
})


setwd("~/MBsnANA/HBana/SCENICana/SCENICOut/")

##load data
load("HUMHB_NB_Neurons.RData")
rn_rna <- row.names(plot.data)
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

##splitting RNA cells
rcells <- chunk2(rn_rna,5)
rlen <- lengths(rcells)

plot.data$Split <- "ND"
rn_rna2 <- rn_rna
for (i in seq(length(rlen))) {
  x <- rlen[i]
  rn_rnax <- sample(rn_rna2, x)
  plot.data[rn_rnax,"Split"] <- letters[i]
  rn_rna2 <- rn_rna2[!(rn_rna2 %in% rn_rnax)]
  
}
table(plot.data$Annotation,plot.data$Split)
save(plot.data, file = "plotdataNBNeusplit.RData")
