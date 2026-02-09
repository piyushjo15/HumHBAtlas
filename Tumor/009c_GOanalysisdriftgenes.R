
suppressPackageStartupMessages({
   library(tidyverse)
    library(Matrix)
  library(ggplot2)
library(gprofiler2)
    library(enrichR)
})
setwd("~/RECORDR/")



##score density PA vs DMG
pa <- read.delim("Scores_PA.txt")
dmg <- read.delim("Scores_DMG.txt")

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
dmg <- read.delim("Scores_DMG.txt",row.names = 1)
dmg <- row.names(dmg)[1:500]

pa <- read.delim("Scores_PA.txt",row.names = 1)
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

                   

