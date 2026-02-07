## merge adjacency matrix from each class run and filter for non-regulation
suppressPackageStartupMessages({
  library(tidyverse)
})

Cls <- commandArgs(trailingOnly = TRUE)

setwd("~/SCENICOut/Class/")
print(paste0("Mergeing adjacency matrix for class: ",Cls))
### ----
adj_a <- read.delim(paste0(Cls,"/cor",Cls,"_A.tsv"))
adj_b <- read.delim(paste0(Cls,"/cor",Cls,"_B.tsv"))

##remove non regulation and low rho
adj_a <- adj_a[!is.na(adj_a$rho),]
adj_a <- adj_a[adj_a$regulation!=0,]

adj_b <- adj_b[!is.na(adj_b$rho),]
adj_b <- adj_b[adj_b$regulation!=0,]

## tfs
tfs_a <- unique(adj_a$TF)
tfs_b <- unique(adj_b$TF)
tfs <- intersect(tfs_a,tfs_b)

##
adj <- c()
for(x in tfs){
  del_a <- adj_a[adj_a$TF==x,]
  del_b <- adj_b[adj_b$TF==x,]
  
  merged_df <- rbind(del_a,del_b) %>%
    group_by(target) %>%
    summarise(rho = mean(rho, na.rm = TRUE),
              importance=mean(importance, na.rm = TRUE),
              count = n(), .groups = "drop")
  ## connection is positive and present in both
  merged_df_fil <- merged_df %>% filter(count > 1, rho > 0.03)
  merged_df_fil$TF <- x
  merged_df_fil$regulation <-1
  adj <- rbind(adj,merged_df_fil)
  rm(merged_df_fil,del_a,del_b,merged_df)
}
rm(x)
adj <- adj[,c("TF","target","importance","regulation","rho")]
write.table(adj, file = paste0(Cls,"/cor",Cls,"_fil.tsv"),
            sep = "\t",quote = FALSE,row.names = FALSE)
print(paste0("Finished fixing adjacency matrix for class: ",Cls))
q()
#