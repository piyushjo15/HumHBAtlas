## this script finds TAD associated with SEs and then identify which genes lies
## in thsoe TAD and assign them to the SEs
suppressPackageStartupMessages({
  library(tidyverse)
})
DIR <- "~/MBsnANA/HBana/ATACana/Outs/Superenhancer/"
setwd(DIR)
Cls <- commandArgs(trailingOnly = TRUE)
print(paste0("Performing analysis for class: ",Cls))
# ## part 1. Find TAD and SE pairs -------------
# ## After obtaining at max top 150 SE per class, I am identifying the
# ## TADs in which these SE lie or overlap
# load(paste0("rose/",Cls,"_SE_cre.RData")) 
# SEs <- data.frame(A=names(se_cre_list))
# SEs <- SEs %>% separate(A,c("chr","start","end"))
# SEs$start <- as.integer(SEs$start)
# SEs$end <- as.integer(SEs$end)
# all_se <- row.names(SEs) <- SEs$Name <- names(se_cre_list)
# #head(SEs)
# ## read appropriate TAD
# ## classes that potentially belongs to Germinal Zone
# cls1 <- c("Progenitors","Glioblast","PreO","GCUBC",
#           "Neuroblasts_A","MES_VAS","Neuroblasts_C")
# if(Cls %in% cls1){
#   tad <- read.delim("GSE77565_GZ_TAD.hg38.sort.bed",header = FALSE)
# }else {
#   tad <- read.delim("GSE77565_CP_TAD.hg38.sort.bed",header = FALSE)
# }
# colnames(tad) <- c("chr","start","end")
# #head(tad)
# 
# ## TO associate TAD(s) associated with SE, I am checking if the
# ## start of SE is greater than start of TAD and end is less than
# ## end of TAD, sometimes an SE can overlap one or more TAD
# tad4se <- list()
# for(se in all_se){
# 
#   chr <- SEs[se,"chr"]
#   st <- SEs[se,"start"]
#   ed <- SEs[se,"end"]
#   tad_sel <- tad[tad$chr==chr,]
#   ## filtering tads that end after SE start and start before SE end
#   tad_sel <- tad_sel[tad_sel$end>st & tad_sel$start<ed,]
#   ## the first TAD is the one that start before SE start
#   if(dim(tad_sel)[1]>0){
#     tad4se[[se]] <- tad_sel
#   }
#   rm(chr,st,ed,tad_sel)
# }
# 
# #save(tad4se, SEs,file = paste0("SEgenes/",Cls,"_tadse.RData"))
# #q()
# 
## part 2. Find TAD, SE and gene combinations ------------
## now after finding the TAD associated with each SE, I find the genes associated
## with each tad and assign them to genes
#load(paste0("SEgenes/",Cls,"_tadse.RData"))

all_se <- names(tad4se)
gene_tss <- read.delim("~/MBsnANA/HBana/ATACana/Outs/extrafiles/gencdH38p13r37CR_genesTSS.bed",header = FALSE)
colnames(gene_tss)[1:4] <- c("chr","start","end","gene")
#head(gene_tss)
genes4se <- list()
for(se in all_se){

  ## TAD(s) associated with SE
  tad <- tad4se[[se]]
  ## TAD boundaries associated with SE
  chr <- SEs[se,"chr"]
  st <- min(tad$start)
  ed <- min(tad$end)

  gene_tss_sel <- gene_tss[gene_tss$chr==chr,]
  gene_tss_sel <- gene_tss_sel[gene_tss_sel$start>st & gene_tss_sel$end <ed,]
  if(dim(gene_tss_sel)[1]>0){
    genes4se[[se]] <- unique(gene_tss_sel$gene)
  }
  rm(tad,chr,st,ed,gene_tss_sel)

}
rm(se)
#table(lengths(genes4se))
save(tad4se, SEs,genes4se, file = paste0("SEgenes/",Cls,"_tadsegenes.RData"))
q()

# ## part 5. Finding SE associated genes based on constituent peaks------------
# ## This is better appraoch
# suppressPackageStartupMessages({
#   library(rGREAT)
#   library(GenomicRanges)
# 
# })
# ## load peaks that constitutes each SE of the class
# 
# load(paste0("rose/",Cls,"_SE_cre.RData"))
# load(paste0("SEgenes/",Cls,"_tadsegenes.RData")) ### genes associated with TAD
# ## convert list to data frame
# all_se <- names(se_cre_list)
# 
# SEcre <- c()
# for(se in all_se){
#   cre <- data.frame(A=se_cre_list[[se]])
#   cre <- cre %>% separate(A,c("chr","start","end"))
#   cre$SE <- se
#   SEcre <- rbind(SEcre,cre)
#   rm(cre)
# }
# rm(se)
# SEcre$strand <- "."
# 
# ###converting to GRanges
# PS_gr <- makeGRangesFromDataFrame(SEcre,
#                          keep.extra.columns=TRUE,
#                          ignore.strand=FALSE,
#                          seqinfo=NULL,
#                          seqnames.field=c("chr"),
#                          start.field="start",
#                          end.field="end",
#                          strand.field="strand",
#                          starts.in.df.are.0based=TRUE)
# ## in order to provide dummy gene-set for the rGREAT to run
# ## selecting some genes and obtianing ENSG id for those genes
# some_genes <- unique(unlist(genes4se))[1:10]
# ens <- read.delim("~/DATA/HB/ann/gencdH38p13r37CR_genes.txt",header = FALSE)
# ens <- ens %>% separate(V1,c("A"),sep="\\.")
# colnames(ens) <- c("Gene.ID","Gene_name")
# ens <- ens[!duplicated(ens$Gene.ID),]
# row.names(ens) <- ens$Gene.ID
# ##dummy geneset to avoid GO analysis
# ens2 <- unique(ens[ens$Gene_name%in%some_genes,"Gene.ID"])
# gs <- list(What=ens2)
# 
# res <- great(PS_gr,gene_sets = gs,
#              tss_source = "Gencode_v37",
#              min_gene_set_size = 1,
#              verbose = FALSE)
# 
# #plotRegionGeneAssociations(res)
# 
# R2G <- getRegionGeneAssociations(res)
# R2G <- data.frame(R2G)
# row.names(R2G) <- paste0(R2G$seqnames,":",R2G$start-1,"-",R2G$end) ##-1 in Peaks
# 
# R2G$Gene.ID <- R2G$Genes <- NA
# SEgreat <- list()
# for(se in all_se){
#   R2G_sel <- R2G[R2G$SE==se,]
#   peaks <- row.names(R2G_sel)
#   gene_sel <- c()
#   for(x in peaks){
#     gx <- unlist(R2G[x,"annotated_genes"])
#     gene_names <- ens[gx,"Gene_name"]
#     gene_sel <- c(gene_sel,gene_names)
#     gene_names <- paste0(gene_names,collapse = ",")
#     gx <- paste0(gx,collapse = ",")
#     R2G[x,"Gene.ID"] <- gx
#     R2G[x,"Genes"] <- gene_names
#     rm(gx, gene_names)
#   }
#   SEgreat[[se]] <- unique(gene_sel)
#   rm(x,peaks,gene_sel)
# }
# 
# #head(R2G)
# save(R2G,SEgreat,genes4se,SEs, file = paste0("SEgenes/",Cls,"_SEgeneGREATv2.RData"))
# q()
## part 6. FInding genes that are identified as target by GREAT and are in TAD in v2 approach-----------
## Now using great analysis and TAD info to find a list of genes associated with
## class specific SEs
load(paste0("SEgenes/",Cls,"_SEgeneGREATv2.RData"))

all_se <- row.names(SEs)
genes <- c()
for(se in all_se){
  tad_genes <- genes4se[[se]]
  great_genes <- SEgreat[[se]]
  if(length(tad_genes)>1){
    del <- intersect(tad_genes,great_genes)
    genes <- c(genes,del)
  }else{
    genes <- c(genes,great_genes)
  }
  rm(tad_genes,great_genes)
}
rm(se)

genes <- unique(genes)
genes2 <- unique(unlist(genes4se))
write(genes, file = paste0("SEgenes/",Cls,"_SEgenev2.txt"))
write(genes2, file = paste0("SEgenes/",Cls,"_SEgeneonlyTADv2.txt"))

q()
########## This another attempt I used to associate genes with SE, upper one is better##############
# ## part3. Use GREAT for SE gene association------
# ## annotate peaks by GREAT 
# suppressPackageStartupMessages({
#   library(rGREAT)
#   library(GenomicRanges)
# 
# })
# ## load data
# load(paste0("SEgenes/",Cls,"_tadsegenes.RData"))
# SEs$strand <- "."
# ###converting to GRanges
# peaks <- row.names(SEs)
# print("Duplicated SE peaks:")
# table(duplicated(peaks))##all unique!!
# PS_gr <- makeGRangesFromDataFrame(SEs,
#                          keep.extra.columns=FALSE,
#                          ignore.strand=FALSE,
#                          seqinfo=NULL,
#                          seqnames.field=c("chr"),
#                          start.field="start",
#                          end.field="end",
#                          strand.field="strand",
#                          starts.in.df.are.0based=TRUE)
# ## in order to provide dummy gene-set for the rGREAT to run
# ## selecting some genes and obtianing ENSG id for those genes
# some_genes <- unique(unlist(genes4se))[1:10]
# ens <- read.delim("~/DATA/HB/ann/gencdH38p13r37CR_genes.txt",header = FALSE)
# ens <- ens %>% separate(V1,c("A"),sep="\\.")
# colnames(ens) <- c("Gene.ID","Gene_name")
# ens <- ens[!duplicated(ens$Gene.ID),]
# row.names(ens) <- ens$Gene.ID
# ##dummy geneset to avoid GO analysis
# ens2 <- unique(ens[ens$Gene_name%in%some_genes,"Gene.ID"])
# gs <- list(What=ens2)
# 
# res <- great(PS_gr,gene_sets = gs,
#              tss_source = "Gencode_v37",
#              min_gene_set_size = 1,
#              verbose = FALSE)
# 
# #plotRegionGeneAssociations(res)
# 
# R2G <- getRegionGeneAssociations(res)
# R2G <- data.frame(R2G)
# row.names(R2G) <- paste0(R2G$seqnames,":",R2G$start-1,"-",R2G$end) ##-1 in Peaks
# #head(R2G)
# print("Peaks annotated by GREAT:")
# 
# table(row.names(R2G) %in% row.names(SEs))
# 
# peaks <- row.names(R2G) ## these are SEs that are associated with GREAT analysis
# R2G$Gene.ID <- R2G$Genes <- NA
# SEgreat <- list()
# for(x in peaks){
#   gx <- unlist(R2G[x,"annotated_genes"])
#   gene_names <- ens[gx,"Gene_name"]
#   SEgreat[[x]] <- gene_names
#   gene_names <- paste0(gene_names,collapse = ",")
#   gx <- paste0(gx,collapse = ",")
#   R2G[x,"Gene.ID"] <- gx
#   R2G[x,"Genes"] <- gene_names
#   rm(gx, gene_names)
# }
# rm(x)
# #head(R2G)
# save(R2G,SEgreat,genes4se,SEs, file = paste0("SEgenes/",Cls,"_SEgeneGREAT.RData"))
# q()
# ## part4. Finding genes that are identified as target by GREAT and are in TAD-----------
# ## Now using great analysis and TAD info to find a list of genes associated with
# ## class specific SEs
# load(paste0("SEgenes/",Cls,"_SEgeneGREAT.RData"))
# 
# all_se <- row.names(SEs)
# genes <- c()
# for(se in all_se){
#   tad_genes <- genes4se[[se]]
#   great_genes <- SEgreat[[se]]
#   if(length(tad_genes)>1){
#     del <- intersect(tad_genes,great_genes)
#     genes <- c(genes,del)
#   }else{
#     genes <- c(genes,great_genes)
#   }
#   rm(tad_genes,great_genes)
# }
# rm(se)
# 
# genes <- unique(genes)
# write(genes, file = paste0("SEgenes/",Cls,"_SEgene.txt"))
# q()
# 



####################### NOT doing this ################
## not doing this --------
# ## This is for annotating the peak region as exonic, intronic or intergenic
# ###ChIPPeakannotation-----------
# 
# dd<- GenomicFeatures::makeTxDbFromGFF("~/DATA/HB/ann/gencdH38p13r37CR.gtf", format = "gtf",
#                      dataSource="Gencodep13r37",
#                      organism = "Homo sapiens")
# 
# 
# out <- ChIPpeakAnno::genomicElementDistribution(PS_gr,
#                                   TxDb = dd,
#                                   #TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#                                   promoterRegion=c(upstream=2000, downstream=500),
#                                   geneDownstream=c(upstream=0, downstream=2000),
#                                   promoterLevel=list(
#                                     # from 5' -> 3', fixed precedence 3' -> 5'
#                                     breaks = c(-2000, -1000, -500, 0, 500),
#                                     labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
#                                                "upstream <500b", "TSS - 500b"),
#                                     colors = c("#FFE5CC", "#FFCA99",
#                                                "#FFAD65", "#FF8E32")))
# out2 <- ChIPseeker::annotatePeak(PS_gr, 
#                                  tssRegion=c(-500, 500),
#                                  TxDb=dd)
# ANN <- data.frame(out$peaks)
# row.names(ANN) <- paste0(ANN$seqnames,":",ANN$start-1,"-",ANN$end)
# table(row.names(ANN) %in% peaks)
# ANN <- out$peaks
# ANN2 <- data.frame(out2@anno)
# row.names(ANN2) <- names(out2@anno)
# save(PS, PS_gr, R2G,ANN,ANN2, file=paste0(DIR_PO,"Peaks_robust_GREAT_annotated.RData")) ##association GREAT
# q()
#### -------------
###########
#