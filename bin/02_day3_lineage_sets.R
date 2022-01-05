# include -----------------------------------------------------------------
library(scater)
library(scran)
library(Seurat)

library(tidyverse)
library(ggthemes)
set.seed(2510)
source('./bin/useful_functions.R')
# input

sce <- readRDS('./data/day_3_processed/ctrl.rds')
sce <- seurat_clusters(sce, ndimred = 20, resolution_param = 1.5) # clusters used to resolve between mk and ery CMPs

# Isolation of high-pseudotime clusters of cells -----------------------------------------------
# Define tip clusters isolating the highest pseudotime scores for cells in each branch
sce$dpt_cluster = '0'

# stem cluster
mask <- sce$mk_gmp_branch2 == 1 
sce$dpt_cluster[mask][order(sce$dpt[mask])][1:300] = '1'

# ly cluster
mask <- sce$mk_gmp_branch2 == 4
sce$dpt_cluster[mask][order(sce$dpt[mask], decreasing = T)][1:300] = '4'

# gmp cluster
mask <- sce$mk_gmp_branch2 == 3 
sce$dpt_cluster[mask][order(sce$dpt[mask], decreasing = T)][1:300] = '3'

# mep cluster
mask <- sce$mk_gmp_branch2 == 2
sce$dpt_cluster[mask][order(sce$dpt[mask], decreasing = T)][1:300] = '2'

# Isolation of mk and erythroid-specific sets
# mk and ery cells are very close. Clusters obtained using the louvain algorithm are used to distinguish between them
tmp <- findMarkers(sce,sce$cluster, direction = 'up', test.type = 'wilcox')

# top markers
tmp2 <- tmp$`12`%>% as_tibble(rownames = 'ID') %>% mutate(Symbol = rowData(sce)[ID,'Symbol'])  %>% 
  dplyr::select(Symbol, FDR, AUC.10) %>% arrange(FDR) # erythroid - basophil signature# megakaryocytic signature
tmp3 <- tmp$`10`%>% as_tibble(rownames = 'ID') %>% mutate(Symbol = rowData(sce)[ID,'Symbol'])  %>%
  dplyr::select(Symbol, FDR, AUC.12) %>% arrange(FDR) # erythroid - basophil signature

# define tip clusters further resolving mk and ery
sce$dpt_cluster2 <- sce$dpt_cluster
# remove dpt_cluster 2 to split it in 5 (ery) and 6 (mk)
sce$dpt_cluster2[sce$dpt_cluster ==2] = '0'

# ery
mask <- sce$cluster == 10
sce$dpt_cluster2[mask][order(sce$dpt[mask], decreasing = T)][1:300] = '5'

# mk
mask <- sce$cluster == 12
sce$dpt_cluster2[mask][order(sce$dpt[mask], decreasing = T)][1:300] = '6'


# lineage markers detection -----------------------------------------------
## part 1: joint mk and ery signature (named mep)
markers <- findMarkers(sce, sce$dpt_cluster, direction = 'up', test.type = 'wilcox', pval.type = 'all', row.data = rowData(sce))
columns = paste0('AUC.', as.character(c(1,2,3,4)))

lin_sets = lapply(c(1,2,3,4), function(x) {
  # cast the x as character 
  name = paste0('AUC.',as.character(x))
  # make up the name of the columns pasting strings
  df =   markers[[as.character(x)]] %>% #+1 due to the clustering starting at 0
    as.data.frame() %>% 
    rownames_to_column(var = 'gene')
  
  df[,name] = 10
  
  df$min = apply(df[,columns],MARGIN = 1,FUN = min)
  df = df[df$FDR < 5e-2, ]
  df[,name] = 0
  df = df[,c('gene','Symbol', 'FDR', columns)]
  return(df)
})
sapply(lin_sets, nrow)  # number of detected markers

# save
names(lin_sets) = c('stem','mep','my','ly')
for (set in 1:length(lin_sets)){
  write.table(format(x = lin_sets[[set]] %>% dplyr::mutate(Symbol = rowData(sce)[gene,'Symbol']), digits =2),
              paste0('./data/day_3_processed/lin_sets/',
                     names(lin_sets)[set],'_all.csv'),
              sep = '\t',quote = F, row.names = F)}

## part 2: separate mk and ery signature (named mep)
markers <- findMarkers(sce, sce$dpt_cluster2, direction = 'up', test.type = 'wilcox', pval.type = 'some', min.prop = 0.75, row.data = rowData(sce))
columns = paste0('AUC.', as.character(c(1,3,4,5,6)))

lin_sets = lapply(c(1,3,4,5,6), function(x) {
  # cast the x as character 
  name = paste0('AUC.',as.character(x))
  # make up the name of the columns pasting strings
  df =   markers[[as.character(x)]] %>% #+1 due to the clustering starting at 0
    as.data.frame() %>% 
    rownames_to_column(var = 'gene')
  
  df[,name] = 10
  
  df$min = apply(df[,columns],MARGIN = 1,FUN = min)
  df = df[df$FDR < 5e-2, ]
  df[,name] = 0
  df = df[,c('gene','Symbol', 'FDR', columns)]
  return(df)
})
sapply(lin_sets, nrow)

# save mk and ery sets
names(lin_sets) = c('stem','my','ly','ery','mk')
for (set in c('ery','mk')){
  write.table(format(x = lin_sets[[set]] %>% dplyr::mutate(Symbol = rowData(sce)[gene,'Symbol']), digits =2),
              paste0('./data/day_3_processed/lin_sets/',
                     set,'_all.csv'),
              sep = '\t',quote = F, row.names = F)}