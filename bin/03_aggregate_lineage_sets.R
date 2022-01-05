library(tidyverse)
library(org.Mm.eg.db)
library(AnnotationHub)

cc.genes <- AnnotationDbi::select(org.Mm.eg.db, keys="GO:0007049", keytype="GOALL", column="ENSEMBL")$ENSEMBL
ribosomial.genes <- AnnotationDbi::select(org.Mm.eg.db, keys="GO:0005840", keytype= "GOALL", column="ENSEMBL")$ENSEMBL

# combine markers sets to obtain lineage sets
lineages <- c('stem', 'ly', 'my', 'mk', 'ery')

top_markers <- lapply(lineages, FUN = function(name){
  d1_markers <- read.table(paste0('./data/day_1_processed/lin_sets/',name,'_all.csv'),sep = '\t', header = T) %>% 
    mutate(rank = rownames(.)) %>% dplyr::select(gene, rank, Symbol) 
  d3_markers <- read.table(paste0('./data/day_3_processed/lin_sets/',name,'_all.csv'),sep = '\t', header = T) %>% 
    mutate(rank = rownames(.)) %>% dplyr::select(gene, rank, Symbol)
  
  # merge
  markers <- full_join(d1_markers, d3_markers, by = 'gene')
  markers$rank.x[is.na(markers$rank.x)] <- nrow(d1_markers) + 1
  markers$rank.y[is.na(markers$rank.y)] <- nrow(d3_markers) + 1
  markers$avg_rank = (as.integer(markers$rank.x) + as.integer(markers$rank.y))/2
  markers <- markers[order(markers$avg_rank), c('gene','Symbol.x')]
  
  # exclude cycle
  markers <- markers[!markers$gene %in% cc.genes,]
  # exclude ribosomial
  markers <- markers[!markers$gene %in% ribosomial.genes,]
  # exclude pseudogenes
  annot <- AnnotationDbi::select(org.Mm.eg.db, keys= as.character(markers$gene), columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
  pseudogenes <- annot$ENSEMBL[grep('pseudogene', annot$GENENAME)]
  markers <- markers[!markers$gene %in% pseudogenes,]

  top_100 <- markers$gene[1:100]
  return(top_100)})

names(top_markers) <- lineages

for(i in lineages){
  line = c(i,  top_markers[[i]],'\n')
  cat(line, file = './data/lin_sets/lineage.sets.gmt', sep ='\t',append = T)}

mep.line <- c('mep', top_markers$mk[1:50], top_markers$ery[1:50], '\n')

cat(mep.line, file = './data/lin_sets/lineage.sets.gmt', sep ='\t',append = T)


##### Gene Symbol 
top_markers <- lapply(lineages, FUN = function(name){
  d1_markers <- read.table(paste0('./data/day_1_processed/lin_sets/',name,'_all.csv'),sep = '\t', header = T) %>% 
    mutate(rank = rownames(.)) %>% dplyr::select(gene, rank, Symbol) 
  d3_markers <- read.table(paste0('./data/day_3_processed/lin_sets/',name,'_all.csv'),sep = '\t', header = T) %>% 
    mutate(rank = rownames(.)) %>% dplyr::select(gene, rank, Symbol)
  
  # merge
  markers <- full_join(d1_markers, d3_markers, by = 'gene')
  markers$rank.x[is.na(markers$rank.x)] <- nrow(d1_markers) + 1
  markers$rank.y[is.na(markers$rank.y)] <- nrow(d3_markers) + 1
  markers$avg_rank = (as.integer(markers$rank.x) + as.integer(markers$rank.y))/2
  markers <- markers[order(markers$avg_rank), c('gene','Symbol.x')]
  
  # exclude cycle
  markers <- markers[!markers$gene %in% cc.genes,]
  # exclude ribosomial
  markers <- markers[!markers$gene %in% ribosomial.genes,]
  
  # exclude pseudogenes
  annot <- AnnotationDbi::select(org.Mm.eg.db, keys= as.character(markers$gene), columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
  pseudogenes <- annot$ENSEMBL[grep('pseudogene', annot$GENENAME)]
  markers <- markers[!markers$gene %in% pseudogenes,]
  
  top_100 <- as.character(markers$Symbol.x[1:100])
  return(top_100)})

names(top_markers) <- lineages

for(i in lineages){
  line = c(i,  top_markers[[i]],'\n')
  cat(line, file = './data/lin_sets/lineage.sets_symbol.gmt', sep ='\t',append = T)}

mep.line <- c('mep', top_markers$mk[1:50], top_markers$ery[1:50], '\n')

cat(mep.line, file = './data/lin_sets/lineage.sets_symbol.gmt', sep ='\t',append = T)



