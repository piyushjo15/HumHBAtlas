suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
    library(ggrepel)
    library(ggplot2)
    library(clusterProfiler)
})
setwd("~/DIPG/")


### PA -----------
sams <- c("DMG01","DMG02","DMG03","DMG04","DMG05")
GSEA <- list()

##astro
astro_up <- list()
astro_down <- list()
for(x in sams){
    load(paste0(x,"/degAstrovsTumAstro_glm_new.RData"))
    keep1 <- res$adj_pval<1e-5
    keep2 <- abs(res$lfc)>=2
    res$significant <- "No"
    res[keep1 & keep2,"significant"] <- "Yes"
    astro_up[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc>0]
     astro_down[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc<0]
    rm(res,keep1,keep2)
}
lengths(astro_up)
lengths(astro_down)

astro_up <- data.frame(table(unlist(astro_up)))
astro_down <- data.frame(table(unlist(astro_down)))
GSEA[["astro_up"]] <- as.character(astro_up$Var1)[astro_up$Freq>2]
GSEA[["astro_down"]] <- as.character(astro_down$Var1)[astro_down$Freq>2]


#oligo
oligo_up <- list()
oligo_down <- list()
for(x in sams){
    load(paste0(x,"/degOPCvsTumOPC_glm_new.RData"))
    keep1 <- res$adj_pval<1e-5
    keep2 <- abs(res$lfc)>=2
    res$significant <- "No"
    res[keep1 & keep2,"significant"] <- "Yes"
    oligo_up[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc>0]
     oligo_down[[x]] <- row.names(res)[res$significant=="Yes" & res$lfc<0]
    rm(res,keep1,keep2)
}
lengths(oligo_up)
lengths(oligo_down)

oligo_up <- data.frame(table(unlist(oligo_up)))
oligo_down <- data.frame(table(unlist(oligo_down)))
GSEA[["oligo_up"]] <- as.character(oligo_up$Var1)[oligo_up$Freq>2]
GSEA[["oligo_down"]] <- as.character(oligo_down$Var1)[oligo_down$Freq>2]


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
save(GSEA, file="DIPG_degs.RData")

### GO analysis
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

save(terms_list,terms_df,file="DIPG_Reac.RData")

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
save(aa,file="DIPG_gprofiler.RData")
q()



