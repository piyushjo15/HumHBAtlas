##This script calculates enrichment scores (using Seurat algorithm) for a gene-set
## for hindbrain data 


suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(tidyverse)
  library(AUCell)
  library(BiocParallel)
})
addmodscore <- function(obj, features, name = "Cluster",
                        pool = row.names(obj), nbin=24,
                        ctrl = 100, k =FALSE, seed = 123){
  # Find how many gene lists were provided. In this case just one.
  cluster.length <- length(x = features)
  # For all genes, get the average expression across all cells (named vector)
  data.avg <- Matrix::rowMeans(x = obj[pool, ])
  # Order genes from lowest average expression to highest average expression
  data.avg <- data.avg[order(data.avg)]
  
  # Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations.
  # The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
  data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                  n = nbin,
                                  labels = FALSE,
                                  right = FALSE)
  # Set the names of the cuts as the gene names
  names(x = data.cut) <- names(x = data.avg)
  # Create an empty list the same length as the number of input gene sets.
  # This will contain the names of the control genes
  ctrl.use <- vector(mode = "list", length = cluster.length)
  
  # For each of the input gene lists:
  for (i in 1:cluster.length) {
    # Get the gene names from the input gene set as a character vector
    features.use <- features[[i]]
    
    # Loop through the provided genes (1:num_genes) and for each gene,
    # find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
    for (j in 1:length(x = features.use)) {
      # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number.
      # We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                          size = ctrl,
                                          replace = FALSE)))
    }
  }
  # Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin,
  # there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling
  # the same gene more than once.
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  
  
  ## Get control gene scores
  
  # Create an empty matrix with dimensions;
  # number of rows equal to the number of gene sets (just one here)
  # number of columns equal to number of cells in input Seurat obj
  ctrl.scores <- matrix(data = numeric(length = 1L),
                        nrow = length(x = ctrl.use),
                        ncol = ncol(x = obj))
  
  # Loop through each provided gene set and add to the empty matrix the mean expression of the control genes in each cell
  for (i in 1:length(ctrl.use)) {
    # Get control gene names as a vector
    features.use <- ctrl.use[[i]]
    # For each cell, calculate the mean expression of *all* of the control genes
    ctrl.scores[i, ] <- Matrix::colMeans(x = obj[features.use,])
  }
  
  ## Get scores for input gene sets
  # Similar to the above, create an empty matrix
  features.scores <- matrix(data = numeric(length = 1L),
                            nrow = cluster.length,
                            ncol = ncol(x = obj))
  
  # Loop through input gene sets and calculate the mean expression of these genes for each cell
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- obj[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  
  # Subtract the control scores from the feature scores -
  # the idea is that if there is no enrichment of the genes in the geneset in a cell,
  # then the result of this subtraction should be ~ 0
  features.scores.use <- features.scores - ctrl.scores
  
  # Name the result the "name" variable + whatever the position the geneset was in the input list, e.g. "Cluster1"
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  
  # Change the matrix from wide to long
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  
  # Give the rows of the matrix, the names of the cells
  rownames(x = features.scores.use) <- colnames(x = obj)
  
  return(features.scores.use)
  
}

##Directories 

##1. load required data -------
### 1.1 load gene-sets ------
## any gene-set, it has to be a list of character arrays
## the name of the list has to be GSEA
## Gene-sets
load("extrafiles/Reactome_Signaling_Pathways.RData")
##TF-GRNs
#load("extrafiles/AllGRNs.RData")
#GSEA <- all_grns
### 1.2 load expression and metadata ------

## this is compiled log transformed expression data for the entire hindbrain
load("lgcountsHB.RData")

##load metadata, combined
load("plotdataHBv2.RData")
## subset by class
# tb <- data.frame(table(plot.data$Class))
# ncells <- min(tb$Freq) ## minimum number of cells are in CBINT
# ncells <- 9000
# ## selecting 9000 cells per class
# cls <- unique(plot.data$Class)
# cells <- c()
# for(x in cls){
#   pddel <- plot.data[plot.data$Class==x,]
#   del_cells <- sample(row.names(pddel),ncells)
#   cells <- c(cells,del_cells)
#   rm(pddel,del_cells)
# }
# rm(x)
# write(cells, file = "CellformodscoreHBv2.txt")

## subset by Cluster
# ## susbet metdata so each class have same numebr of cells
# #tb <- data.frame(table(plot.data$Cluster))
# ncells <- 200 ## taking at max cells per cluster
# ## selecting 200 cells per cluster
# cls <- unique(plot.data$Cluster)
# cells <- c()
# for(x in cls){
#   del_cells <- row.names(plot.data)[plot.data$Cluster==x]
#   if(length(del_cells)>ncells){
#     del_cells <- sample(del_cells,ncells)
#   }
#   cells <- c(cells,del_cells)
#   rm(del_cells)
# }
# rm(x)
# write(cells, file = "CellformodscoreHBv3.txt")

cells <- readLines("CellformodscoreHBv2.txt")

plot.data <- plot.data[row.names(plot.data) %in% cells,]

mdt <- mdt[,row.names(plot.data)]

# ### 1.3 remove ribosomal and mitochndiral genes from genesets--------
#load("~/MBsnANA/Diemana/scepro/genesnoMTRibo.RData")
#GSEAx <- GSEA
#GSEA <- list()
#cls <- names(GSEAx)
#for(x in cls){
#  cand <- GSEAx[[x]]
#  cand <- cand[cand %in% genes]
#  #to remove small sets, if needed
#  if(length(cand)>19){
#    GSEA[[x]] <- cand
#  }
#  rm(cand)
#}

## 2. module score----
scr <- addmodscore(mdt,features = GSEA)
colnames(scr) <- names(GSEA)
row.names(scr) <- row.names(plot.data)

save(scr, plot.data, file = "AUC_MS/pdHB_GeneSet_MS.RData")

rm(scr)
q()
####################### Analysis #############
## New AUC bubble plot 
## Using enrichment odd ratio and proportion of a cells beloning to a 
## class or cluster in top 90/95 percentile
## The odd ratio is how many cells for a class are above a cut-off (let say 5%)
## compared to the same for others
## I will also down sample for this

load("AUC_MS/pdHB_GeneSet_MS.RData")

## for class
cl_lv <- c("Progenitors","Neuroblasts_I","Neuroblasts_II","Neuroblasts_III",
            "GCUBC","GC","CBINT","Purkinje", "Neurons_I", "Neurons_II","Glioblasts",
            "PreO","MOG","MixedGlia","MES_VAS","Immune")
scr[1:4,1:4]
mps <- colnames(scr) ### this is the names of sets
auc_rank_mat <- apply(scr, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC
table(row.names(plot.data)==row.names(auc_rank_mat))

# Set high-AUC threshold (global or per-signature)
auc_co <- 0.90 # or 0.95 when comparing all clusters in the Hb 

# Create result list
res_list <- list()

for (x in mps) {
  auc_vec <- auc_rank_mat[, x]
  high_vec <- auc_vec > auc_co
  
  ## create a datafrom
  df <- data.frame(cell = rownames(auc_rank_mat), 
                   high_auc = high_vec, 
                   cluster = plot.data$Class)
  row.names(df) <- df$cell
  # Total cell above rank cut-off across clusters
  all_high <- sum(df$high_auc)
  ncells <- nrow(df)
  
  res_sig <- df %>%
    group_by(cluster) %>%
    summarise(
      tot = n(), ### total number of cells in a cluster
      a = sum(high_auc), # In cluster, high AUC
      b = tot - a,# In cluster, not high AUC
      c = all_high-a,# Outside cluster, high AUC
      d = (ncells- tot)-c,
      prop = (a / tot)*100,# % of above cut-off cells
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      odds_ratio = (a * d) / (b * c + 1e-10),
      p = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value,
      MP = x
    )
  
  res_list[[x]] <- res_sig
  rm(auc_vec,high_vec,ncells,all_high,res_sig)
}

# Combine and adjust p-values
merge_res <- bind_rows(res_list) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "BH"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-10))

## cu-off here depends on the maximum value of log2OR
merge_res <- merge_res %>%
  group_by(MP) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 5), 0)) %>%
  ungroup()

head(merge_res)
merge_res <- data.frame(merge_res)

merge_res$cluster <- factor(merge_res$cluster,levels = cl_lv)
merge_res$MP <- factor(merge_res$MP,levels = rev(mps))

## for SE 
merge_res <- merge_res %>%
  group_by(MP) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 4), 0)) %>%
  ungroup()


merge_res$cluster <- factor(merge_res$cluster,levels = cl_lv) ## cluster order
merge_res$MP <- factor(merge_res$MP,levels = mps) ## set order

p <- ggplot(merge_res, aes(y = MP, x = cluster)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 10)) +
  scale_color_gradient2(low = "white",mid = "grey70", high = "black",midpoint = 2.5,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.x = element_text(colour = "black",size=14,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(colour = "black",size=14),
        panel.border=element_rect(colour = 'black',linewidth=1),
        axis.title=element_blank())+
scale_y_discrete(position = "right") 



pdf("Figures/LDA18MP_HBlv1_MSenrprop.pdf", width = 8, height = 6, pointsize = 20)
print(p)
dev.off()
