### find upregulated and down regulated genes for PA samples
library(Matrix)

setwd("~/MBsnANA/Diemana/scepro/RECORDR")
sam <- commandArgs(trailing=TRUE)
print(paste0("extract DEG for sample :", sam))
GSEA<- list()

##Astro
load(paste0(sam,"/degAstrovsTumAstro_glm_new.RData"))

res <- data.frame(res)
res$Sig <- ifelse(res$adj_pval < 0.001,"Yes","No")

sig <- list()

sig[["up"]] <- row.names(res)[res$lfc>=2 & res$Sig=="Yes"]
sig[["down"]]  <- row.names(res)[res$lfc<= -2 & res$Sig=="Yes"]

GSEA[["Astro"]] <- sig

rm(res,sig)

##OPC
load(paste0(sam,"/degOPCvsTumOPC_glm_new.RData"))

res <- data.frame(res)
res$Sig <- ifelse(res$adj_pval < 0.001,"Yes","No")

sig <- list()

sig[["up"]] <- row.names(res)[res$lfc>=2 & res$Sig=="Yes"]
sig[["down"]]  <- row.names(res)[res$lfc<= -2 & res$Sig=="Yes"]

GSEA[["OPC"]] <- sig

rm(res,sig)

##Astro
load(paste0(sam,"/degMicrovsTumMicro_glm_new.RData"))

res <- data.frame(res)
res$Sig <- ifelse(res$adj_pval < 0.001,"Yes","No")

sig <- list()

sig[["up"]] <- row.names(res)[res$lfc>=2 & res$Sig=="Yes"]
sig[["down"]]  <- row.names(res)[res$lfc<= -2 & res$Sig=="Yes"]

GSEA[["Micro"]] <- sig

rm(res,sig)
save(GSEA, file=paste0(sam,"/GSEA_deg.RData"))
q()
## now extracting DEG per compnonent across samples
load("PA74/GSEA_deg.RData")
GSEA_PA74 <- GSEA

load("PA87/GSEA_deg.RData")
GSEA_PA87 <- GSEA

load("PA102/GSEA_deg.RData")
GSEA_PA102 <- GSEA
rm(GSEA)

opc_up <- c(GSEA_PA74[["OPC"]][["up"]],GSEA_PA87[["OPC"]][["up"]],GSEA_PA102[["OPC"]][["up"]])
opc_up <- data.frame(table(opc_up))
opc_down <- c(GSEA_PA74[["OPC"]][["down"]],GSEA_PA87[["OPC"]][["down"]],GSEA_PA102[["OPC"]][["down"]])
opc_down <- data.frame(table(opc_down))


astro_up <- c(GSEA_PA74[["Astro"]][["up"]],GSEA_PA87[["Astro"]][["up"]],GSEA_PA102[["Astro"]][["up"]])
astro_up <- data.frame(table(astro_up))
astro_down <- c(GSEA_PA74[["Astro"]][["down"]],GSEA_PA87[["Astro"]][["down"]],GSEA_PA102[["Astro"]][["down"]])
astro_down <- data.frame(table(astro_down))

micro_up <- c(GSEA_PA74[["Micro"]][["up"]],GSEA_PA87[["Micro"]][["up"]],GSEA_PA102[["Micro"]][["up"]])
micro_up <- data.frame(table(micro_up))
micro_down <- c(GSEA_PA74[["Micro"]][["down"]],GSEA_PA87[["Micro"]][["down"]],GSEA_PA102[["Micro"]][["down"]])
micro_down <- data.frame(table(micro_down))

GSEA <- list()
GSEA[["OPC_up"]] <- as.character(opc_up$opc_up[opc_up$Freq>1])
GSEA[["OPC_down"]] <- as.character(opc_down$opc_down[opc_down$Freq>1])
GSEA[["Astro_up"]] <- as.character(astro_up$astro_up[astro_up$Freq>1])
GSEA[["Astro_down"]] <- as.character(astro_down$astro_down[astro_down$Freq>1])
GSEA[["Micro_up"]] <- as.character(micro_up$micro_up[micro_up$Freq>1])
GSEA[["Micro_down"]] <- as.character(micro_down$micro_down[micro_down$Freq>1])
lengths(GSEA)
save(GSEA, file="~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/PA_DEG_GSEA.RData")











