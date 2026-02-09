##merging drft scores from each run of tumor sample with all three normal combinations
## here shown for PA, similar analysis done for DMG
suppressPackageStartupMessages({
  library(tidyverse)
  library(gprofiler2)
})
setwd("~/RECORDR/")


##amir
#PA01
dfa <- read.csv("PA01/gene_context_drift_scores_nor_PA01_a.csv",row.names = 1)
dfb <- read.csv("PA01/gene_context_drift_scores_nor_PA01_b.csv",row.names = 1)
dfc <- read.csv("PA01/gene_context_drift_scores_nor_PA01_c.csv",row.names = 1)
dfa$ID <- row.names(dfa)
dfb$ID <- row.names(dfb)
dfc$ID <- row.names(dfc)

df1 <- list(dfa, dfb, dfc) %>%
  reduce(full_join, by = "ID") %>%
  transmute(
    ID,
    total_score = rowMeans(select(., starts_with("total_score")), na.rm = TRUE)
  )

rm(dfa,dfb,dfc)
row.names(df1) <- df1$ID

#PA02
dfa <- read.csv("PA02/gene_context_drift_scores_nor_PA02_a.csv",row.names = 1)
dfb <- read.csv("PA02/gene_context_drift_scores_nor_PA02_b.csv",row.names = 1)
dfc <- read.csv("PA02/gene_context_drift_scores_nor_PA02_c.csv",row.names = 1)
dfa$ID <- row.names(dfa)
dfb$ID <- row.names(dfb)
dfc$ID <- row.names(dfc)

df2 <- list(dfa, dfb, dfc) %>%
  reduce(full_join, by = "ID") %>%
  transmute(
    ID,
    total_score = rowMeans(select(., starts_with("total_score")), na.rm = TRUE)
  )

rm(dfa,dfb,dfc)
row.names(df2) <- df2$ID

#PA03
dfa <- read.csv("PA03/gene_context_drift_scores_nor_PA03_a.csv",row.names = 1)
dfb <- read.csv("PA03/gene_context_drift_scores_nor_PA03_b.csv",row.names = 1)
dfc <- read.csv("PA03/gene_context_drift_scores_nor_PA03_c.csv",row.names = 1)
dfa$ID <- row.names(dfa)
dfb$ID <- row.names(dfb)
dfc$ID <- row.names(dfc)

df3 <- list(dfa, dfb, dfc) %>%
  reduce(full_join, by = "ID") %>%
  transmute(
    ID,
    total_score = rowMeans(select(., starts_with("total_score")), na.rm = TRUE)
  )
rm(dfa,dfb,dfc)
row.names(df3) <- df3$ID

genes <- c(row.names(df1),row.names(df2),row.names(df3))
genes <- data.frame(table(genes))
genes <- as.character(genes$genes)[genes$Freq>=2]
head(df1)

merged <- data.frame(PA01=df1[genes,"total_score"],PA02=df2[genes,"total_score"],PA03=df3[genes,"total_score"])
row.names(merged) <- genes
merged$Avg <- rowMeans(merged,na.rm = TRUE)
merged <- merged[order(merged$Avg,decreasing = TRUE),]
merged$Rank <- seq(dim(merged)[1])
ggplot(merged,aes(x=Rank,y=Avg))+
  geom_point(size=0.5)+
  geom_hline(yintercept =1,color="red",linetype="dashed")+
  theme_bw()
write.table(merged, file = "Scores_PA.txt",
            sep = "\t",col.names = TRUE, row.names = TRUE,quote = FALSE)

