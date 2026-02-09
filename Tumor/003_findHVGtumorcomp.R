### findHVG in PA, PFA and DMG for DEG
suppressPackageStartupMessages({
  library(scran)
  library(scater)
 
})
## first finding for PA----------
DIR_OUT="~/PA"
setwd(DIR_OUT)
load("lgcountsPA.RData")
load("plotdataPA_ANN.RData")

## find HVG for OPC
pd <- plot.data[plot.data$ANN2 %in% c("OPC","OPC_Astro"),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)

##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "PA_OPC_HVG.RData")
rm(pd,mdtx,dec,dec1,dec2,hvg_tum)

## find HVG for Astro
pd <- plot.data[plot.data$ANN2 %in% c("Astro"),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "PA_Astro_HVG.RData")
rm(pd,mdtx,dec,dec1,dec2,hvg_tum)

## find HVG for Microg
pd <- plot.data[grep("Micro",plot.data$ANN3),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "PA_Micro_HVG.RData")

rm(pd,mdtx,dec,dec1,dec2,hvg_tum)
rm(plot.data,mdt)


## then finding for DMG----------
DIR_OUT="~/DIPG/"
setwd(DIR_OUT)
load("lgcountsDIPG.RData")
load("plotdataDIPG_ANN.RData")



## find HVG for OPC
pd <- plot.data[plot.data$ANN2 %in% c("Oligo"),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "DMG_OPC_HVG.RData")
rm(pd,mdtx,dec,dec1,dec2,hvg_tum)

## find HVG for Astro
pd <- plot.data[plot.data$ANN2 %in% c("Astro"),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "DMG_Astro_HVG.RData")
rm(pd,mdtx,dec,dec1,dec2,hvg_tum)

## find HVG for Microg
pd <- plot.data[grep("Micro",plot.data$ANN3),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "DMG_Micro_HVG.RData")

rm(pd,mdtx,dec,dec1,dec2,hvg_tum)
rm(plot.data,mdt)


## then finding for DMG----------
DIR_OUT="~/PFA/"
setwd(DIR_OUT)
load("lgcountsGJPFA.RData")
load("plotdataGJPFA_ANN.RData")

load("~/extrafiles/genesnoMTRibo.RData")
genes1 <- genes
genes <- row.names(mdt)
genes <- intersect(genes1,genes)
## find HVG for Ependymal
pd <- plot.data[plot.data$ANN2 %in% c("Ependymal"),]

mdtx <- mdt[genes,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "PFA_Epen_HVG.RData")
rm(pd,mdtx,dec,dec1,dec2,hvg_tum)

## find HVG for Microg
pd <- plot.data[grep("Micro",plot.data$ANN3),]
mdtx <- mdt[,row.names(pd)]
dec<- modelGeneVar(mdtx)
##remvoe unintersting genes
dec1 <- dec[!is.na(dec$p.value),]
dec2 <- dec[dec$mean>0.02,]
hvg_tum <- getTopHVGs(dec2,n=3000)
hvg_tum <- hvg_tum[!is.na(hvg_tum)]

save(hvg_tum, file = "PFA_Micro_HVG.RData")

rm(pd,mdtx,dec,dec1,dec2,hvg_tum)
rm(plot.data,mdt)

q()