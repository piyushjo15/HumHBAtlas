suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
    library(ggrepel)
    library(ggplot2)
    library(clusterProfiler)
})
setwd("~/MBsnANA/Diemana/scepro/RECORDR/")


### PA -----------
sams <- c("MUV013","MUV014","MUV051")
GSEA <- list()

##astro
epen_up <- list()
epen_down <- list()
for(x in sams){
    load(paste0(x,"/degEpenvsTumEpen_glm_new.RData"))
    keep1 <- res$adj_pval<1e-5
    keep2 <- abs(res$lfc)>=2
    res$significant <- "No"
    res[keep1 & keep2,"significant"] <- "Yes"
    epen_up[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc>0]
     epen_down[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc<0]
    rm(res,keep1,keep2)
}
lengths(epen_up)
lengths(epen_down)

epen_up <- data.frame(table(unlist(epen_up)))
epen_down <- data.frame(table(unlist(epen_down)))
GSEA[["epen_up"]] <- as.character(epen_up$Var1)[epen_up$Freq>1]
GSEA[["epen_down"]] <- as.character(epen_down$Var1)[epen_down$Freq>1]

#micro
micro_up <- list()
micro_down <- list()
for(x in sams){
    load(paste0(x,"/degMicrovsTumMicro_glm_new.RData"))
    keep1 <- res$adj_pval<1e-5
    keep2 <- abs(res$lfc)>=2
    res$significant <- "No"
    res[keep1 & keep2,"significant"] <- "Yes"
    micro_up[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc>0]
     micro_down[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc<0]
    rm(res,keep1,keep2)
}
lengths(micro_up)
lengths(micro_down)

micro_up <- data.frame(table(unlist(micro_up)))
micro_down <- data.frame(table(unlist(micro_down)))
GSEA[["micro_up"]] <- as.character(micro_up$Var1)[micro_up$Freq>2]
GSEA[["micro_down"]] <- as.character(micro_down$Var1)[micro_down$Freq>2]
##save
lengths(GSEA)
save(GSEA, file="GojoPFA_degs.RData")

### GO analysis
#load("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/PID_KEGG_REAC_WIKI.RData")
load("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/PID_GOBP_REAC_WIKI.RData")

msigdbr_t2g = pid_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
#msigdbr_t2g <-msigdbr_t2g[grep("REAC",msigdbr_t2g$gs_name),]
head(msigdbr_t2g)
clones <- names(GSEA)
terms_list <- c()
terms_df <- c()

for(x in clones){
  ks <- enricher(gene = GSEA[[x]], TERM2GENE = msigdbr_t2g)
  res <- ks@result
  keep <- res$p.adjust<=0.05 ##0.05 for wiki, 0.01 for react
  if(any(keep)){
    res <- res[keep,c("ID","GeneRatio","geneID","p.adjust")]
    row.names(res) <- NULL
    res$Cluster <- x
    #terms_list[[x]] <- res$ID
    terms_df <- rbind(terms_df,res)
    rm(ks,keep,res)
  }
  
}

save(terms_list,terms_df,file="GojoPFA_Reac.RData")

library(gprofiler2)
gostres <- gost(query =GSEA, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO:BP","REAC"), as_short_link = FALSE, highlight = FALSE)
aa <- gostres$result
aa <- aa[,c("query","term_id","p_value","source","term_name","parents")]
save(aa,file="GojoPFA_gprofiler.RData")
q()



