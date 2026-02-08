library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")

## Obtain key names from MergeMotifs
setwd("~/DeepHB/Outs/")
pat_names <- readLines("all_pat_names.txt")
pat_names <- data.frame(X=pat_names)
pat_names <- pat_names %>% separate(X,c("A","B"))
pat_names <- pat_names$B

repl_python()
import pickle
with open('pm_lstmv2_57MP_73_.pkl', 'rb') as f:
  pattern_matrix=pickle.load( f)

exit
pattern_matrix<- py$pattern_matrix
colnames(pattern_matrix) <- paste0("seqid_",sprintf("%03d", 1:165))
tops <- readLines("TopicMPs.txt")
row.names(pattern_matrix) <- tops
pattern_matrix[1:4,1:4]
save(pattern_matrix, file = "~/DeepHB/Outs/Seqletlogcounts_allPat.RData")
