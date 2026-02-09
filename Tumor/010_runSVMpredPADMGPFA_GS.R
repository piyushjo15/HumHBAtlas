## SVM classfier for tumor cells using normal hindbrain as a reference 
## using a Grid search for best parameters
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")

suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(batchelor)
})
DIR_INP <- "~/RNA/"
setwd(DIR_INP)

##load and prepare reference data -----
load("plotdataHBv2.RData")

plot.data_ref <- plot.data


load("lgcountsHB.RData")
mdt_ref <- mdt[,row.names(plot.data_ref)]
genes <- row.names(mdt)
rm(mdt)
genes1 <- row.names(mdt_ref)
dec <- modelGeneVar(mdt_ref)
hvg_ref <- getTopHVGs(dec)
load("~/extrafiles/genesnoRPMT.RData")
hvg_ref <- hvg_ref[hvg_ref %in% genes]
load("~/extrafiles/HGNCsexchr.RData")
hvg_ref <- hvg_ref[!(hvg_ref %in% SX)]
hvg_ref <- hvg_ref[1:4000]
Cluster <- plot.data_ref$Cluster
mdt_ref <- cosineNorm(mdt_ref[hvg_ref,])

##load and prepare test data -----

##PA
load("~/PA/lgcountsPA.RData")
mdt_pa <- mdt

load("~/PA/topHVGPAnoRPMTSX.RData")
hvg_pa <- com.hvg
load("~/PA/plotdataPA_ANN.RData")
pd_pa <- plot.data[,c("Batch","ANN1","ANN2","ANN3","ind_cluster")]
pd_pa$Tumour <- "PA"
mdt_pa <- mdt_pa[,row.names(pd_pa)]

## DMG
load("~/DIPG/lgcountsDIPG.RData")
mdt_dmg <- mdt
load("~/DIPG/topHVGDIPGnoRPMTSX.RData")
hvg_dmg <- com.hvg
load("~/DIPG/plotdataDIPG_ANN.RData")
pd_dmg <- plot.data[,c("Batch","ANN1","ANN2","ANN3","ind_cluster")]
pd_dmg$Tumour <- "DMG"
mdt_dmg <- mdt_dmg[,row.names(pd_dmg)]

## PFA
load("~/PFA/lgcountsGJPFA.RData")
mdt_gj <- mdt
genes2 <- row.names(mdt_gj)

load("~/PFA/topHVGGJPFAnoRPMTSX_1500.RData")
hvg_gj <- com.hvg

load("~/PFA/plotdataGJPFA_ANN.RData")
pd_gj <- plot.data[,c("Batch","ANN1","ANN2","ANN3","ind_cluster")]
pd_gj$Tumour <- "PFA_GJ"

mdt_gj <- mdt_gj[,row.names(pd_gj)]

## combine targets
genes  <- intersect(genes,genes2)
hvg_tar <- unique(c(hvg_pa,hvg_dmg,hvg_gj))
hvg_tar <- intersect(hvg_tar,genes)

mdt_pa <- cosineNorm(mdt_pa[hvg_tar,])
mdt_dmg <- cosineNorm(mdt_dmg[hvg_tar,])
mdt_gj <- cosineNorm(mdt_gj[hvg_tar,])

mdt_tar <- cbind(mdt_pa,mdt_dmg,mdt_gj)
plot.data <- rbind(pd_pa,pd_dmg,pd_gj)

hvg <- intersect(hvg_ref,hvg_tar)
paste0("the number of genes used in SVM model are:")
print(length(hvg))
mdt_ref <- mdt_ref[hvg,]
mdt_tar <- mdt_tar[hvg,]

## conver genes to ENSG id
ens <- read.delim("~/extrafiles/ENSG2HGNC.txt")
hvg_ref2 <- ens[hvg,"GeneID"]
write(hvg_ref2, file = "~/extrafiles/reference_genes_glial.txt")
row.names(mdt_ref) <- row.names(mdt_tar) <- hvg_ref2

mdt_ref <- t(mdt_ref)
mdt_tar <- t(mdt_tar)


##run grid search on SVM
repl_python()

import numpy as np
from scipy import sparse
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import make_scorer, balanced_accuracy_score
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import joblib

X = r.mdt_ref
if not sparse.isspmatrix_csr(X):
  X = X.tocsr()

X = X.astype(np.float32)

X_test = r.mdt_tar
if sparse.isspmatrix(X_test):
  X_test = X_test.tocsr().astype(np.float32)
else:
  X_test = np.asarray(X_test, dtype=np.float32)

y = np.asarray(r.Cluster)

svc = LinearSVC(penalty="l2", loss="squared_hinge",
                     class_weight="balanced", dual=False,
                     max_iter=20000, random_state=0)

param_grid = {"C": [0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

# Balanced accuracy is a good choice for imbalanced multiclass;
# you could use "accuracy" if classes are fairly balanced
search = GridSearchCV(
  svc,
  param_grid=param_grid,
  scoring="balanced_accuracy",
  cv=cv, n_jobs=4, refit=True, verbose=1
)
search.fit(X, y)

best_params_out = search.best_params_
best_cv_balacc = float(search.best_score_)

# --- calibrate with CV (no leakage) ---
# rebuild an unfitted pipeline with best C
bestC = best_params_out["C"]
base = LinearSVC(penalty="l2", loss="squared_hinge",
                 class_weight="balanced", dual=False,
                 C=bestC, max_iter=20000, random_state=0)

cal = CalibratedClassifierCV(base, method="isotonic", cv=cv)  
cal.fit(X, y)

##saving calibrated model
joblib.dump({
  "model": cal,
  "best_params": best_params_out,
  "best_cv_balacc": best_cv_balacc
}, "/home/extrafiles//HB_linear_svc_calibrated_glial.pkl")

# --- predict ---
test_l1  = cal.predict(X_test)
test_prob = cal.predict_proba(X_test)

print("Best params (LinearSVC):", best_params_out)
print("Best CV balanced_accuracy:", best_cv_balacc)

exit

a <- py$test_l1
b <- py$test_prob

plot.data$SVM_prilab_ANN2 <- plot.data$SVM_finlab_ANN2 <- a
prob <- apply(b, 1, max)
plot.data$SVM_prob_ANN2 <- prob
plot.data$SVM_finlab_ANN2[prob < 0.5] <- "ND"


save(a,b,plot.data, file="SVMpred_PADMGPFA_HBref.RData")
q()

load("SVMpred_WNTSHH_HBref.RData")
d <- length(unique(plot.data$SVM_prilab_ANN2))
ggplot(plot.data,aes(x=Batch,fill=SVM_prilab_ANN2))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = colorRampPalette(pal_igv()(51))(d))+
  theme_minimal()
  theme(legend.position = "None")




