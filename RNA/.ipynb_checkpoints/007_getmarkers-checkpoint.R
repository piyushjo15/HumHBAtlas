# find differentially expressed genes per class

suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
})

load("lgcountsHB.RData")
load("plotdataHBv2.RData")
## subset
mdt <- mdt[,row.names(plot.data)]
##using marker for all
markers_w <- findMarkers(mdt, groups = plot.data$Class, test.type="wilcox",
                         pval.type="some")
markers_b <- findMarkers(mdt, groups = plot.data$Class, test.type="binom",
                         pval.type="some")
markers_t <- findMarkers(mdt, groups = plot.data$Class, test.type="t",
                         pval.type="some")

save(markers_w,markers_b,markers_t, file = "markersClassHB.RData")
cls <- names(markers_w)
GSEA_w <- list()
GSEA_t <- list()
for(x in cls){
    del <- markers_w[[x]]
    GSEA_w[[x]] <- row.names(del)[del$summary.AUC>0.65 & del$FDR<1e-5]
    del <- markers_t[[x]]
    GSEA_t[[x]] <- row.names(del)[del$summary.logFC>.5 & del$FDR<1e-3]
    rm(del)
}
lengths(GSEA_w)
lengths(GSEA_t)

save(GSEA_w,GSEA_t, file="markerGSEAclass.RData")
q()
