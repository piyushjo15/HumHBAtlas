library(universalmotif)
library(ggseqlogo)
DIR <- "~/DeepHB/"
setwd(DIR)

jaspar <- read_homer("~/DeepHB/Outs/allpat_lstmv2_57MP_73.motif")

motifs_names <- unlist(lapply(jaspar, function(m){return(m@name)}))
motifs_len <- unlist(lapply(jaspar, function(m){return(length(unlist(str_split(m@consensus,""))))}))
motifs_consensus <- unlist(lapply(jaspar, function(m){return(m@consensus)}))
df <- data.frame(Name=motifs_names,Con=motifs_consensus)

write.table(df, file="~/DeepHB/Outs/allpat_lstmv2_57MP_735.txt",
            sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE)

for(i in 1:length(motifs_names)){
  w=motifs_len[i]*.4
  pdf(paste0("DeepHBlogos/",motifs_names[i],".pdf"),width = w,height = 2,pointsize = 1)
  print(ggseqlogo(jaspar[[i]]@motif, method = "bits"))
  dev.off()
}


ggseqlogo(jaspar[[i]]@motif, method = "bits")
