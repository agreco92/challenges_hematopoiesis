# this script contains the analysis of Tie2 subsets response to 
# sepsis in terms of lineage markers, proliferation

set.seed(2607)

library(Matrix)
library(scran)
library(scater)
library(Seurat)
library(org.Mm.eg.db)

library(tidyverse)
library(patchwork)
library(ggthemes)

sce = readRDS('../data/day_1_processed/sce2.rds')
sce = sce[,sce$population == 'HSC']

rownames(sce) = uniquifyFeatureNames(ID = rowData(sce)[['ID']], names = rowData(sce)[['Symbol']])

# re-normalize data using pooling method
qckCl = quickCluster(sce, min.size = 100, use.ranks =T, method = 'igraph', d=10 )
sce = computeSumFactors(sce, cluster = qckCl) # remove nothing
logcounts(sce) <- normalizeCounts(sce, size_factors = sizeFactors(sce))


# Response across Tie2 subsets --------------------------------------------------------------------------------

# annotate Tie2 subsets
sce$Tie2_p = counts(sce)['Tek',] > 0

# differential expression analysis for each subset
mks_Tie2p = findMarkers(sce[,sce$Tie2_p & sce$population=='HSC'],
                        groups = sce$condition[sce$Tie2_p & sce$population=='HSC'])

mks_Tie2n = findMarkers(sce[,(!sce$Tie2_p) & sce$population=='HSC'],
                        groups = sce$condition[(!sce$Tie2_p) & sce$population=='HSC'])

df = full_join((mks_Tie2p$sepsis %>% as.data.frame() %>% rownames_to_column('gene') %>% select(gene, FDR, logFC.ctrl)),
               (mks_Tie2n$sepsis %>% as.data.frame() %>% rownames_to_column('gene') %>% select(gene, FDR, logFC.ctrl)), 
               by = 'gene', suffix = c('.tek_p','.tek_n'))


# load lineage markers
load('../../lineage_sets/list_geneSets_tier2.RData')
x = x[!grepl('immat_ery', names(x))] # remove gene set associated with ery-baso split

# plot logFC for Tie2 subsets
df$Symbol = rowData(sce)[df$gene,'Symbol']
df$lin_marker = ''

for(lin in names(x)){
  df$lin_marker[df$Symbol %in% x[[lin]]] = lin
}

p1 <- ggplot(df %>% filter(FDR.tek_p < 0.01 |FDR.tek_n < 0.01),aes(x = logFC.ctrl.tek_p, y = logFC.ctrl.tek_n)) + 
  geom_point(col = 'grey60', size =1, stroke =0, alpha =.3) +
  geom_abline(slope =1, intercept = 0) + 
  geom_point(data = df %>% filter(lin_marker != ''), aes(col = lin_marker), size =1.75, stroke=0) + 
  scale_color_tableau() + theme_bw() + 
  labs(title  = 'rho ~ 0.94',x = 'logFC in Tie2+ subset', y = 'logFC in Tie2- subset')

p1

# compute correlation
df %>% filter(FDR.tek_p < 0.01 |FDR.tek_n < 0.01) %>% summarise(rho = cor(logFC.ctrl.tek_p,logFC.ctrl.tek_n))

# cell cycle analysis -----------------------------------------------------

# import cell cycle genes
s.genes <- readRDS('../../docs/seurat_cyclegenes_mouse/s.genes')
g2m.genes <- readRDS('../../docs/seurat_cyclegenes_mouse/g2m.genes')

gc()

# import in Seurat object
seu = as.Seurat(sce)

# score cell cycle
seu = CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)

## cell cycle embedding ----------------------------------------------------

# embed cells using only cell cycle genes

seu$sum_log = log2(seu$sum)
seu_cc = seu[c(s.genes,g2m.genes), ]

seu_cc = ScaleData(seu_cc)
seu_cc@reductions = list()
seu_cc = RunPCA(seu_cc, features = rownames(seu_cc))
seu_cc = RunUMAP(seu_cc, reduction = 'pca', dims = 1:5)

seu_cc@active.ident = seu_cc$Phase %>% as.factor()

seu@reductions$umap_cc = CreateDimReducObject(embeddings = seu_cc@reductions$umap@cell.embeddings, key = 'umap_cc')

p <- FeaturePlot(seu_cc, features = c('S.Score','G2M.Score')) & theme(axis.text = element_blank(), legend.position = 'none')
p

##  assignment of phase using principal graphs --------------------------------------------------------

### fit graph to embedding ----------------------------------------------------------------------


library("ElPiGraph.R")
library(RColorBrewer)
X = seu_cc@reductions$umap@cell.embeddings
circle_epg <- computeElasticPrincipalCircle(X = X, NumNodes = 30, Lambda = 1e-05, Mu =.1, Do_PCA = F)

p1 <- PlotPG(X = X, TargetPG = circle_epg[[1]],
             NodeLabels = 1:nrow(circle_epg[[1]]$NodePositions),
             LabMult = 2.5, PointSize = NA, p.alpha = .1)[[1]] + scale_y_reverse() + theme_bw() + 
  theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = 'none') + 
  labs(title = 'graph fitting')
p1

circle_graph = ConstructGraph(circle_epg[[1]])
circle_e2e <- GetSubGraph(Net = circle_graph, Structure = 'circle', Circular = T, Nodes = 30)

PartStruct <- PartitionData(X = X, NodePositions = circle_epg[[1]]$NodePositions)
ProjStruct <- project_point_onto_graph(X = X,
                                       NodePositions = circle_epg[[1]]$NodePositions,
                                       Edges = circle_epg[[1]]$Edges$Edges,
                                       Partition = PartStruct$Partition)

### annotate nodes ------------------------------------------------------------------------------

df = data.frame(node = PartStruct$Partition)
df[,c('S.Score','G2M.Score')] = seu@meta.data[,c('S.Score','G2M.Score')]
df[,c('lib.size')] = seu@meta.data[,c('sum_log')]

p2 = df %>% group_by(node) %>% summarise(across(everything(), mean)) %>% 
  column_to_rownames('node') %>% scale() %>%
  pheatmap::pheatmap(color = rev(colorRampPalette(colors = brewer.pal(11,'Spectral'))(101)),
                     cutree_rows = 12, cluster_cols = F, treeheight_row = 20)
#pheatmap::pheatmap()
p2 <- p2[[4]]

p <- wrap_plots(p1,p2, widths = c(0.75,0.5))
p # used to annotate

cc_phase = plyr::mapvalues(df$node, from = seq(1,30), 
                           to = c('G1','G0/G1','S','post_mito_G1','S','G2M','G0/G1','post_mito_G1','G1','G0/G1',
                                  'G1','G2M','G0/G1','G0/G1','G1','G0/G1','G2M','post_mito_G1', 'G1','G0/G1',
                                  'G0/G1','G2M','G1','G0/G1','post_mito_G1','G0/G1','S','G0/G1','G0/G1','G0/G1'))
seu_cc$cc_phase = factor(cc_phase, levels = c('G0/G1', 'G1','S','G2M','post_mito_G1'))

seu_cc@active.ident = seu_cc$cc_phase

as.data.frame(seu_cc@reductions$umap@cell.embeddings) %>% mutate(phase = seu_cc$cc_phase, G2M = seu_cc$G2M.Score, S = seu_cc$S.Score)
ggplot(df, aes(x = G2M, y = S)) + geom_point(aes(col = phase))

## figure 
df = data.frame(cc_phase = factor(cc_phase, levels = c('G0/G1', 'G1','S','G2M','post_mito_G1')), 
                Tie2 = seu$Tie2_p, condition = seu$condition)
df$Tie2 = plyr::mapvalues(df$Tie2, from = c(TRUE, FALSE), to = c('Tie2+', 'Tie2-'))
df$cycling = df$cc_phase != 'G0/G1'
df$condition = gsub('ctrl','control', df$condition)

p <-  df %>% mutate(cycling = cc_phase != 'G0/G1') %>% 
  group_by(Tie2, condition) %>% summarize(cycle = mean(cycling)) %>%
  ggplot(aes(x = condition, y = cycle)) + geom_col(aes(fill = Tie2), position = 'dodge') + 
  labs(y = 'cycling portion', fill = 'LT-HSC \n subset') + 
  ylim(c(0,1))

