
# suppressPackageStartupMessages({
#   library(ArchR)
#   library(scran)
#   library(scater)
#   library(BSgenome.Hsapiens.UCSC.hg38)
#   library(monaLisa)
#   library(universalmotif)
#   library(TFBSTools)
#   library(chromVAR)
#   library(chromVARmotifs)
# })
# 
# set.seed(456)
# addArchRThreads(threads = 8)
# 
# #setwd
# #Next, I load the reference genome
# addArchRGenome("hg38")
# DIR <- "~/MBsnANA/HBana/ATACana/Outs/"
# DIRS <- "~/MBsnANA/HBana/ATACana/HBATAC_Scripts/"
# 
# setwd(DIR)
# #
# 
# 
# 
# proj <- loadArchRProject(path = "comATACx/")
# 
# #add a set of background peaks used in computing deviations
# #proj <- addBgdPeaks(proj, force = T)
# 
# ## JASPAR cluster annotaiton
# #pfms_ret <- homerToPFMatrixList("~/MBsnANA/HBana/ATACana/Outs/extrafiles/JASPAR2024_root.motif", n = 100L)
# pfms_ret <- homerToPFMatrixList("~/MBsnANA/HBana/ATACana/Outs/extrafiles/selected_allpat_lstmv2_57MP_73_220725.motif", n = 100L)
# pwms <- toPWM(pfms_ret)
# 
# #jaspar <- read_homer("~/MBsnANA/HBana/ATACana/Outs/extrafiles/JASPAR2024_root.motif")
# #motifs_names <- unlist(lapply(jaspar, function(m){return(m@name)}))
# 
# deephb <- read_homer("~/MBsnANA/HBana/ATACana/Outs/extrafiles/selected_allpat_lstmv2_57MP_73_220725.motif")
# motifs_names <- unlist(lapply(deephb, function(m){return(m@name)}))
# 
# names(pwms) <- motifs_names
# 
# # proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2024_rrot", 
# #                             annoName = "Jaspar", force = T, motifPWMs = pwms)
# proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "DeepHB", 
#                             annoName = "DeepHB", force = T, motifPWMs = pwms)
# 
# proj <- addDeviationsMatrix(
#   ArchRProj = proj, 
#   #peakAnnotation = "Jaspar",
#   peakAnnotation = "DeepHB",
#   matrixName = "DeepHB_MotifMatrix", 
#   force = TRUE,
#   threads=1
# )
# saveArchRProject(ArchRProj = proj)
## saving the peak matrix
# Get Matrix

suppressPackageStartupMessages({
  library(ArchR)
 
})
set.seed(456)
addArchRThreads(threads = 8)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/MBsnANA/HBana/ATACana/Outs/"
setwd(DIR)

proj <- loadArchRProject(path = "comATACx/")
se <- getMatrixFromProject(proj, useMatrix = "DeepHB_MotifMatrix")

proj <- addImputeWeights(proj, reducedDims="IterativeLSI_int",nRep = 4)
iW <- getImputeWeights(proj)

# 
# mat_motif <- imputeMatrix(
#   mat = assay(se),
#   imputeWeights = iW,
#   threads = getArchRThreads(),
#   verbose = FALSE,
#   logFile = createLogFile("imputeMatrix")
# )
# save(mat_motif, file = "extrafiles/Imputed_motif_mat_selpat.RData")


mat_motif <- imputeMatrix(
  mat = assay(se,"z"),
  imputeWeights = iW,
  threads = getArchRThreads(),
  verbose = FALSE,
  logFile = createLogFile("imputeMatrix")
)
save(mat_motif, file = "extrafiles/Imputed_motif_matZ_selpat.RData")
q()

