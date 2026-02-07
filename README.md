# HumHBAtlas
This repository contains code related to developing human hindbrain multi-omics single-cell atlas.

Following folders contains scripts for specific analyses:
1. RNA: Pre- and post-processing single-nucleus RNA-seq data, integrating entire transcriptomic atlas, clustering, and identify meta-gene programs.
2. ATAC: Pre- and post-processing single-nucelus ATAC-seq data, intergrating entire chromatin accessbility data, label trasnfer from transcriptomic atlas, peak calling and identify meta-regulatory program.
3. DeepHB: Deeplearning based syntax indetification from meta-regulatory program
4. SCENIC: Integrating transcriptomic and chromatin acceessbility profiles per class to obtain TF-GRNs across classes, perform downstream analysis and identify context dependent TF activity across classes.
5. Pontine: TF-GRN analysis comparing pontine nuclei lineage to granule cell lineage.
6. Tumor: Pre- and post-processing tumor data, performing label trasnfer using hindbrain transcriptomic atlas, perform diffential gene expression analysis and perform context drift analysis
7. Misc: Miscellanous scrips
8. extrafiles: Some relevant files used in analysis.
