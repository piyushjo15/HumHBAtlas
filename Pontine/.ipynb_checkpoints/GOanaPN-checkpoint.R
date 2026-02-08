## performing GO analysis for gene-set associated with Pontine Neurons lineage
suppressPackageStartupMessages({
  library(enrichR)
  library(tidyverse)
  library(msigdbr)
  library(clusterProfiler)
})

DIR_GRN <- "~/MBsnANA/HBana/SCENICana/SCENICOut/PN/"

## 

load(paste0(DIR_GRN,"/metadata/PN_gmt_v2.RData"))
GSEA <- regs
### MSigDB
pid_gene_sets = msigdbr(species = "Homo sapiens", category = "C2" #positional geneset
pid_gene_sets = msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")

head(pid_gene_sets)
msigdbr_t2g = pid_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

cls <- names(GSEA)
terms <- c()
for(x in cls){
  ks <- enricher(gene = GSEA[[x]], TERM2GENE = msigdbr_t2g)
  res <- ks@result
  keep <- res$p.adjust<=0.05
  if(any(keep)){
    res <- res[keep,c("ID","GeneRatio","p.adjust","geneID")]
    row.names(res) <- NULL
    res$Cluster <- x
    terms <- rbind(terms,res)
    rm(ks,keep,res)
  }
  
}
rm(x)

