# includes ----------------------------------------------------------------
library(scater)
library(scran)
library(Seurat)

library(tidyverse)
library(ggthemes)
set.seed(2510)
source('./bin/useful_functions.R')

sce <- readRDS('./data/day_3_processed/sce.rds')
sce <- sce[,sce$condition == 'sepsis']

# dimensionality reduction ------------------------------------------------
sce <- runPCA(sce, subset_row =  sce@int_metadata$hvg_no_cc)
dmap <- destiny::DiffusionMap(reducedDim(sce,'PCA')[,1:20])
reducedDim(sce, 'DiffusionMap') <- dmap@eigenvectors

#
rownames(sce) <- uniquifyFeatureNames(ID = rowData(sce)[,'ID'], names = rowData(sce)[,'Symbol'])

# diffusion pseudotime ----------------------------------------------------
## Isolate tip points
# stem tip
stem_tip = as.integer(which.max(as.logical(sce$population == 'HSC') * reducedDim(sce,'DiffusionMap')[,1]))
sce$stem_tip = FALSE
sce$stem_tip[stem_tip] = TRUE

# mk tip
mk_tip = as.integer(which.max(as.logical(sce$population == 'CMP') * reducedDim(sce,'DiffusionMap')[,4]))
sce$mk_tip = FALSE
sce$mk_tip[mk_tip] = T

# ery tip
ery_tip = as.integer(which.max(as.logical(sce$population == 'CMP') * reducedDim(sce,'DiffusionMap')[,7]))
sce$ery_tip =FALSE
sce$ery_tip[ery_tip] = T

# ly tip
ly_tip = as.integer(which.min(as.logical(sce$population == 'MPP') * reducedDim(sce,'DiffusionMap')[,3]))
sce$ly_tip =FALSE
sce$ly_tip[ly_tip] = T

# my tip
my_tip = as.integer(which.max(as.logical(sce$population == 'CMP') * reducedDim(sce,'DiffusionMap')[,3]))
sce$my_tip = FALSE
sce$my_tip[my_tip] = T

## dpt, branching
dpt <- destiny::DPT(dmap, tips = c(stem_tip, mk_tip, my_tip))
plot(dpt, col_by = 'branch')
# substitute NA values with -1
sce$mk_gmp_branch = -1
tmp <- dpt@branch[,1]
tmp[is.na(tmp)] <- -1
sce$mk_gmp_branch = tmp
sce$mk_gmp_branch <- as.character(sce$mk_gmp_branch)

# pseudotime from stem tip cell
sce$dpt = dpt$DPT4574

# save Single Cell Experiment object
saveRDS(sce, './data/day_3_processed/sepsis.rds')