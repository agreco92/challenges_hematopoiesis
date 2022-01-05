# include -----------------------------------------------------------------
library(tidyverse)
library(ggthemes)
library(scater)
library(scran)
library(DropletUtils)
library(AnnotationHub)
library(org.Mm.eg.db)
set.seed(2510)
# preprocessing -----------------------------------------------------------
# read 10x output
folders <- paste0(dir('./data/day_1',full.names = T),'/10x/')
sce <- read10xCounts(folders)
sce$Sample <- gsub('./data/day_1/|/10x/','',sce$Sample)
sce$population <- factor(gsub('.*_','',sce$Sample), levels = c('HSC','ST','MPP','CMP'))
sce$condition <- 'ctrl'
sce$condition[grep('sepsis',sce$Sample)] <- 'sepsis'

# quality control - filter out low-quality cells
## Identify mitochondrial transcripts
ens.mm <- AnnotationHub()[["AH64461"]]
location <- mapIds(ens.mm, keys=rownames(sce), keytype="GENEID", column="SEQNAME")
is.mito <- which(location=="MT")

## Display low-quality cells per sample
tmp <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))
reasons <- quickPerCellQC(tmp, percent_subsets=c("subsets_Mito_percent"), batch = sce$Sample)
colSums(as.matrix(reasons))
reasons$sample <- sce$Sample
reasons %>% as_tibble() %>% group_by(sample) %>% summarise(sum(discard))

colData(sce)[,colnames(tmp)] <- tmp
sce$discard <- reasons$discard

gridExtra::grid.arrange(
  plotColData(sce, x="Sample", y="sum", colour_by="discard", point_size =.5) +  scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="Sample", y="detected", colour_by="discard" ,point_size =.5)  + scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="Sample", y="subsets_Mito_percent", colour_by="discard", point_size =.5)  + ggtitle("Mito percent"),
  ncol=1)

## discard low-quality cells
sce <- sce[,!reasons$discard]

# normalization -----------------------------------------------------------
# scran pooling method with pre-clustering
tmp <- quickCluster(sce, min.size = 500)
sce <- computeSumFactors(sce, cluster = tmp)
sce <- logNormCounts(sce)

# feature selection -------------------------------------------------------
genevar <- modelGeneVar(sce, block = sce$condition)
tmp <- metadata(genevar$per.block$ctrl)
hvg <- getTopHVGs(genevar, var.threshold = 1e-2)

# exclude cell cycle genes to exclude cell cycle signal
cc.genes <- AnnotationDbi::select(org.Mm.eg.db, keys="GO:0007049", keytype="GOALL", column="ENSEMBL")$ENSEMBL
length(intersect(cc.genes, hvg))
hvg_cc_out <- hvg[!hvg%in%cc.genes]

sce@int_metadata$hvg <- hvg
sce@int_metadata$hvg_no_cc <- hvg_cc_out

# dimensionality reduction ------------------------------------------------
sce <- runPCA(sce, subset_row = hvg_cc_out)

# elbow plot
percent.var <- attr(reducedDim(sce), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")

dmap <- destiny::DiffusionMap(reducedDim(sce,'PCA'))
reducedDim(sce, 'DiffusionMap') <- dmap@eigenvectors

# remove outliers in Diffusion coordinates
# WARNING: despite specifying the seed, the diffusion map will produce slightly different (but equivalent) outcomes 
# (see ??scater::calculateDiffusionMap), thus requiring to customize code 
# to eliminate outliers based on the outcome obtained at the time of execution
plotDiffusionMap(sce, col = 'population', ncomponents = 3, shape = 'population') # diagnostic plot to identify outliers
sce$diff_outlier <- as.logical(isOutlier(reducedDim(sce,'DiffusionMap')[,1],type = 'lower',nmads = 5)) | as.logical(isOutlier(reducedDim(sce,'DiffusionMap')[,3],type = 'higher',nmads = 15))
plotDiffusionMap(sce, col = 'diff_outlier', ncomponents = c(3,2), shape = 'population')
table(sce$diff_outlier) # number of outliers
sce <- sce[,!sce$diff_outlier]

# recompute PCA, diffusion map
sce <- runPCA(sce, subset_row = hvg_cc_out) #recompute PCA
dmap <- destiny::DiffusionMap(reducedDim(sce,'PCA'))
reducedDim(sce, 'DiffusionMap') <- dmap@eigenvectors

# save output -------------------------------------------------------------
saveRDS(sce,'./data/day_1_processed/sce.rds')