## this scripts will get control peaks, 1: highly acessible and 2: random
DIRO <- "~/MBsnANA/HBana/ISM/Outs/TopicMP_CREs/"
setwd(DIRO)
##random background peaks----
rand_peaks <- read.delim("~/MBsnANA/HBana/ISM/Outs/extrafiles/HG38_random_peaks.bed",header = FALSE)
head(rand_peaks)
len <- dim(rand_peaks)[1]
rand_peaks2 <- rand_peaks[sample(seq(len),2000),c(1:3)]
rand_peaks2$V4 <- paste0(rand_peaks2$V1,":",rand_peaks2$V2,"-",rand_peaks2$V3)
rand_peaks2$V5 <- "MP0"
head(rand_peaks2)

write.table(rand_peaks2,file = paste0(DIRO,"MP51_CREs.bed"),
            sep = "\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
##highy accessbile peaks background peaks----
