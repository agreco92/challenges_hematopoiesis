# include -----------------------------------------------------------------
library(scater)
library(scran)
library(Seurat)

library(tidyverse)
library(ggthemes)
set.seed(2510)
source('./bin/useful_functions.R')
# input
sce <- readRDS('./data/day_1_processed/sce.rds')
sce <- sce[,sce$condition == 'ctrl']

# dimensionality reduction ------------------------------------------------
sce <- runPCA(sce, subset_row =  sce@int_metadata$hvg_no_cc)
dmap <- destiny::DiffusionMap(reducedDim(sce,'PCA')[,1:20])
reducedDim(sce, 'DiffusionMap') <- dmap@eigenvectors


# diffusion pseudotime ----------------------------------------------------
## Isolate tip points
# stem tip
stem_tip = as.integer(which.min(as.logical(sce$population == 'HSC') * reducedDim(sce,'DiffusionMap')[,1]))
sce$stem_tip = FALSE
sce$stem_tip[stem_tip] = TRUE

# mk tip
mk_tip = as.integer(which.max(as.logical(sce$population == 'CMP') * reducedDim(sce,'DiffusionMap')[,2]))
sce$mk_tip = FALSE
sce$mk_tip[mk_tip] = T

# ery tip
ery_tip = as.integer(which.min(as.logical(sce$population == 'CMP') * reducedDim(sce,'DiffusionMap')[,4]))
sce$ery_tip =FALSE
sce$ery_tip[ery_tip] = T

# ly tip
ly_tip = as.integer(which.min(as.logical(sce$population == 'MPP') * reducedDim(sce,'DiffusionMap')[,3]))
sce$ly_tip =FALSE
sce$ly_tip[ly_tip] = T

# my tip
my_tip = as.integer(which.min(as.logical(sce$population == 'CMP') * reducedDim(sce,'DiffusionMap')[,2]))
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

# assign lymphoid branch by specifying lymphoid tip point 
dpt2 <- destiny::DPT(dmap, tips = c(stem_tip, ery_tip, ly_tip))
plot(dpt2, col_by = 'branch')
sce$mk_gmp_branch2 <- sce$mk_gmp_branch
sce$mk_gmp_branch2[dpt2@branch[,1]==3] = '4'

# pseudotime from stem tip cell
sce$dpt = dpt$DPT9597

# save Single Cell Experiment object
saveRDS(sce, './data/day_1_processed/ctrl.rds')
