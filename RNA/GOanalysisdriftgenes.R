library(enrichR)

suppressPackageStartupMessages({
   library(tidyverse)
    library(Matrix)
  library(ggplot2)
library(gprofiler2)
})
setwd("~/MBsnANA/Diemana/scepro/RECORDR/")

setwd("~/MBsnANA/Diemana/scepro/RECORDR/")


##score density PA vs DMG
pa <- read.delim("Scores_PA_amir_pr_090126.txt")
dmg <- read.delim("Scores_DMG_amir_pr_090126.txt")

rn <- intersect(row.names(pa),row.names(dmg))
pa <- pa[rn,]
dmg <- dmg[rn,]
df <- data.frame(PA=pa$Avg,DMG=dmg$Avg)

p <- ggplot(df2,aes(x=value,fill=variable))+
  geom_density(alpha=0.8,color = "black")+
  scale_fill_manual(values = c("white","pink3"))+
    theme_bw()+
  theme(axis.text= element_text(size=14),
        axis.title = element_blank(),
        axis.ticks.x = element_line(size=1),
        legend.position = "None",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth =2))+
geom_vline(xintercept = 1,linetype="dashed",linewidth=1)
p
pdf("Figures/PAvsDMG_density.pdf",width = 3,height = 3,pointsize = 10)
print(p)
dev.off()
##gprofielr
dmg <- read.delim("Scores_DMG_amir_pr_090126.txt",row.names = 1)
dmg <- row.names(dmg)[1:500]

pa <- read.delim("Scores_PA_amir_pr.txt",row.names = 1)
pa <- row.names(pa)[1:500]
GSEA <- list()
GSEA[["DMG"]] <- dmg
GSEA[["PA"]] <- pa
rm(dmg,pa)
##gprofiler
gostres <- gost(query =GSEA, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("REAC"), as_short_link = FALSE, highlight = FALSE)
aa <- gostres$result
aa <- aa[,c("query","term_id","p_value","source","intersection_size","term_name","parents")]
head(aa)

clones <- c("PA","DMG")
aa_fil <- c()
GSEA_r <- list()

for(x in clones){
  aad <- aa[aa$query==x,]
  ## lets remove ancestors
  terms <- aad$term_id
  termx <- all_par <-c()
  ## remove a term if it is a parent of any high ranked term
  for(i in 1:length(terms)){
    if(!(terms[i] %in% all_par)){
      termx <- c(termx,terms[i])
    }
    pars <- unlist(aad[i,"parents"])
    all_par <- unique(c(all_par,pars))
    rm(pars)
  }
  rm(i,all_par)
  aady <- aad[aad$term_id %in% termx,]
  rm(termx,terms)
  ## now remove chldren term of a highly ranked term
  terms <- aady$term_id
  termx <- termsy <- c()
  for(i in 1:length(terms)){
    pars <- unlist(aady[i,"parents"])
    keep <- pars %in% termsy
    if( !(any(keep))){
      termx <- c(termx,terms[i])
    }
    termsy <- c(termsy,terms[i])
    rm(pars,keep)
  }
  rm(i)
  aadz <- aady[aady$term_id %in% termx,]
  aa_fil <- rbind(aa_fil,aadz[1:10,])
  
  GSEA_r[[x]] <- aadz$term_id
  rm(aad,aady,terms,termx,termsy,aadz)
  
}

rm(x)
lengths(GSEA_r)
aa_fil$Rank <- rev(seq(10))
aa_fil$term_namex <- paste0(aa_fil$term_name,"_",aa_fil$query)
aa_fil$term_namex <- factor(aa_fil$term_namex,level=rev(aa_fil$term_namex))

p <-  ggplot(aa_fil, aes(x=Rank,y=term_namex, size=intersection_size))+
  geom_point(aes(color=p_value))+
  scale_size("intersection_size", range = c(4,10)) +
  scale_color_viridis_c(direction = -1,option = "C" , 
                        name = "FDR")+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size=14),
        panel.grid.major = element_line(linewidth = 1),
        panel.grid.minor = element_blank())+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25))

                   

ks <- enrichr(gene = genes,databases = c("Reactome_Pathways_2024","WikiPathways_2024_Human","GO_Biological_Process_2025"))

res <- bind_rows(ks)
res2 <- res %>% separate(Overlap,c("Count","B"))
res2 <- res2[res2$Adjusted.P.value<0.05,c("Term","Count","Adjusted.P.value","Genes")]
res2 <- res2[res2$Adjusted.P.value<0.1,c("Term","Count","Adjusted.P.value","Genes")] #PA

res2 <- res2[order(res2$Adjusted.P.value,decreasing = FALSE),]

write.table(res2, file = "DMG_GO.txt",sep = "\t",quote = FALSE)
write.table(res2, file = "PA_GO.txt",sep = "\t",quote = FALSE)

res2 <- read.delim("DMG_GO_fil.txt")
res2 <- read.delim("PA_GO_fil.txt")

res2$Term <- factor(res2$Term,levels = rev(res2$Term))
res2$Count <- as.numeric(as.character(res2$Count))  # handles character/factor
res2$Rank <- dim(res2)[1]+1-seq(dim(res2)[1])

head(res2)
p <- ggplot(res2, aes(x=Rank,y=Term, size=Count))+
  geom_point(aes(color=Adjusted.P.value))+
  scale_size("Count", range = c(4,10)) +
  scale_color_viridis_c(direction = -1,option = "C" , 
                        name = "FDR")+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.ticks = element_line(linewidth = 2),
        axis.text = element_text(size=14),
        panel.grid.major = element_line(linewidth = 2),
        panel.grid.minor = element_blank())+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25))
pdf(paste0("Figures/DMG_ReacWP.pdf"),width = 8,height =7,pointsize = 10)
print(p)
dev.off()
pdf(paste0("Figures/PA_ReacWP.pdf"),width =8,height =7,pointsize = 10)
print(p)
dev.off()


suppressPackageStartupMessages({
  library(enrichR)
  library(tidyverse)
  library(msigdbr)
  library(clusterProfiler)
})


### MSigDB
pid_gene_sets1 = msigdbr(species = "Homo sapiens", category = "C2")#positional geneset
pid_gene_sets2 = msigdbr(species = "Homo sapiens", category = "C5")#positional geneset
pid_gene_sets <- rbind(pid_gene_sets1,pid_gene_sets2)

pid_gene_sets <- pid_gene_sets[pid_gene_sets$gs_subcollection %in% c("CP:WIKIPATHWAYS","CP:REACTOME","GO:BP"),]
head(pid_gene_sets)
msigdbr_t2g = pid_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

ks <- enricher(gene = genes, TERM2GENE = msigdbr_t2g)
res <- ks@result
keep <- res$p.adjust<=0.05

res <- res[keep,]

res <- res[,c("ID","Count","p.adjust","geneID")]
res <- res[order(res$p.adjust,decreasing = FALSE),]

write.table(res2, file = "DMG_GO.txt",sep = "\t",quote = FALSE)
write.table(res2, file = "PA_GO.txt",sep = "\t",quote = FALSE)

res2 <- read.delim("DMG_GO.txt")
res2 <- read.delim("PA_GO.txt")

k <- res$ID
k <- gsub("REACTOME_","",k)
k <- gsub("WP_","",k)
res$ID <- str_to_sentence(k, locale = "en")
res$ID <- factor(res$ID,levels = rev(res$ID))
res$Count <- as.numeric(as.character(res$Count))  # handles character/factor
res$Rank <- dim(res)[1]+1-seq(dim(res)[1])
row.names(res) <- NULL
head(res)
p <- ggplot(res, aes(x=Rank,y=ID, size=Count))+
  geom_point(aes(color=p.adjust))+
  scale_size("Count", range = c(4,10)) +
  scale_color_viridis_c(direction = -1,option = "C" , 
                        name = "FDR")+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.ticks = element_line(linewidth = 2),
        axis.text = element_text(size=14),
        panel.grid.major = element_line(linewidth = 2),
        panel.grid.minor = element_blank())+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25))
pdf(paste0("Figures/DMG_ReacWP.pdf"),width = 8,height =7,pointsize = 10)
print(p)
dev.off()
pdf(paste0("Figures/PA_ReacWP.pdf"),width =8,height =7,pointsize = 10)
print(p)
dev.off()


if(any(keep)){
  res <- res[,c("ID","GeneRatio","p.adjust","geneID")]
  row.names(res) <- NULL
  res$Cluster <- x
  terms <- rbind(terms,res)
  rm(ks,keep,res)
  
}


