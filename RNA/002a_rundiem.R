## DIEM  based true cell calling for each single nuclei RNA seq
## The input is intronic and exnoic alignments obtained from 001_runSTARSolo.sh
# 1. Data input----
suppressPackageStartupMessages({
  library(diem)
  library(ggplot2)
  library(gridExtra)
  library(Matrix)
  library(scater)
  library(scran)

})

#argument
args = commandArgs(trailingOnly=TRUE) ## library id

print(paste0("processing sample id: ",args))
#Load the premrna counts
DIR="rawmatrix/"
counts.in <- read_10x(paste0(DIR,args[1],"/GeneFull"))
dim(counts.in)
# load the exon part

counts.ex <- read_10x(paste0(DIR,args[1],"/Gene"))
dim(counts.ex)
#table(row.names(counts.in)==row.names(counts.ex))

## intronic read fraction plot------
tot <- colSums(counts.in)
exonic <- colSums(counts.ex)
#reads not present in exonic reads are intronic?
intronic <- tot-exonic  
frac <- intronic/tot
mdfr <- data.frame(cbind(frac,tot))

#filtering Nan
mdfr <- mdfr[!(is.na(frac)),]
#sorting by tot reads
mdfr <- mdfr[order(-mdfr$tot),]
#rank
mdfr$Rank <- seq(1:dim(mdfr)[1])
##subsetting 20k cells
mdfr2 <- mdfr[1:20000,]

#find similar cells----

rm(mdfr, intronic, exonic, tot, frac)
#check colnames and rownames
if(all(colnames(counts.in) ==colnames(counts.ex)) & all(row.names(counts.in) == row.names(counts.ex))){
  
  #head(counts.ex)[,1:10]
  genes <- read.delim("extrafiles/gencdH38p13r37CR_genesN.txt", header = FALSE) ##ENS id
  genes <- unique(genes$V1)
  
  #finding genes that were counted less with intronic marix
  mean.in <- rowMeans(counts.in)
  mean.ex <- rowMeans(counts.ex)
  d <- (mean.ex > mean.in)
  print(table(d))
  #replacing above genes' data in intronic matrix from exonic matrix
  counts.in2 <- counts.in[!d,]
  dim(counts.in2)
  counts.ex2 <- counts.ex[d,]
  dim(counts.ex2)
  
  counts.in2 <- rbind(counts.in2,counts.ex2)
  dim(counts.in2)
  counts.in2 <- counts.in2[genes,]
  rm(genes)
}

## 1.1 Creating an SCE object------

mdt <- create_SCE(counts.in2, name=args[1])
dim(mdt)
#class(mdt)
### remoe all count matrices
rm(counts.ex2, counts.ex, counts.in)
#processing using diem

## 1.2 Plotting quality metrics
mt_genes <- grep(pattern = "^MT-", x = rownames(mdt@gene_data),
                 ignore.case = TRUE, value = TRUE)
mdt <- get_gene_pct(x = mdt, genes = mt_genes, name = "pct.mt")
genes <- grep(pattern = "^MALAT1$", x = rownames(mdt@gene_data),
              ignore.case = TRUE, value = TRUE)
mdt <- get_gene_pct(x = mdt, genes = genes, name = "MALAT1")


# head(drop_data)
# summary(drop_data)
#
# 2. Running DIEM-----

# The following outlines the steps involved in `diem`:
#
#   1. Speciyfing the test set and the debris set.
# 2. Running PCA on test set
# 2. Running k-means on the PCs to initialize the clusters.
# 4. Running a EM to estimate the parameters of the Dirichlet-multinomial
# mixture model.
# 5. Classifying droplets based on their likelihood.

## 2.1 Specifying test and debris droplets


mdt <- set_debris_test_set(mdt, top_n = 15000,
                           min_counts = 200,
                           min_genes = 200) ##this will be filtered later

length(mdt@test_set)
length(mdt@bg_set)

## 2.2 Running PCA on the test set and intialization-----

mdt <- filter_genes(mdt, cpm_thresh = 10)
genes <- gene_data(mdt)
#summary(genes)
##default pca parameters, genes=2000, pcs=30
mdt <- get_pcs(mdt, threads = 4)


## 3.3 Initializing clusters
#`k_init = 30` says to initialize k-means with 30 centers,
#`min_size_init = 20` says to take clusters with less than 20 droplets
#and assign them to the next closest cluster.
#`nstart_init = 50`  specify to run k-means 50 times and pick the best run.
#`iter.max_init` says each run can have at most 100 iterations.
# We run multiple starts because k-means may converge to a local optimum.
# The run with the lowest total within sum of
# squares is selected.
mdt <- init(mdt,
            k_init = 30,
            nstart_init = 50,
            min_size_init = 50,
            seedn = 123,
            threads = 4)

## 3.4 Running EM using DM like before-----
mdt <- run_em(mdt, threads=4)
mdt <- assign_clusters(mdt)
#less max_genes make debris score lower, use default settings
mdt <- estimate_dbr_score(mdt)

#de_genes <- debris_genes(mdt)
#head(de_genes)
mdt <- call_targets(mdt,thresh_score = .5, min_genes = 500)
mdt.pd <- droplet_data(mdt)
colnames(mdt.pd)[1:2] <- c("nUMIs", "nGenes")
mdt.pd <- mdt.pd[order(-mdt.pd$nUMIs),]
mdt.pd$Rank <- seq(dim(mdt.pd)[1])
mdt.pd$nUMIs <- mdt.pd$nUMIs+1

# save ----
save(mdt, mdfr2,mdt.pd, file=paste0("diemv4",args[1],".RData"))
q()
