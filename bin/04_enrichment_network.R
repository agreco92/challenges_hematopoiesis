# this script is intended to produce a unified representation of GO terms and  NetPath pathways enriched in the
# overlap between DEGs in HSC,ST,MPP CMP_my

# includes ----------------------------------------------------------------
library(tidyverse)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(igraph)

# import gene sets --------------------------------------------------------
de_list <- as.character(read.csv('./docs/DEGs/hsc_st_mpp_cmp.my_shared_d1_up.csv', header = F)[,1])
go_terms <- readRDS('./data/lin_sets/GO_terms/Mm.c5.bp.v7.1.entrez.rds') # obtained from http://bioinf.wehi.edu.au/MSigDB/v7.1/ because of https://www.biostars.org/p/294026/

netpath <- fgsea::gmtPathways('./data/lin_sets/geneset_db.gmt')
netpath <- netpath[grepl('NetPath', names(netpath))]

# convert to entrez id
de_list <- mapIds(org.Mm.eg.db, keys=de_list, column="ENTREZID", keytype="SYMBOL")

# define gene universe
universe <- unique(c(as.character(unlist(go_terms)), 
                     as.character(de_list),
                     as.character(unlist(netpath))))
universe = paste0('ez_',universe)

# matrix gene by pathway (for each gene, the pathways that contain it)
df <- matrix(F, nrow = length(universe), ncol = 1+ length(go_terms) + length(netpath), 
             dimnames = list(universe, c('de_list', names(go_terms), names(netpath)))) %>%
  as.data.frame()

# fill gene-by-pathway matrix

## DE
df[paste0('ez_',de_list),'de_list'] <- T
## GO
for(n in names(go_terms)){
  df[paste0('ez_',go_terms[[n]]), n] <- T
}
## NETPATH
for(n in names(netpath)){
  df[paste0('ez_',netpath[[n]]), n] <- T
}

# compute enrichment ------------------------------------------------------
# compute p values using Fisher test
p.vals <- rep(1, ncol(df))
names(p.vals) <- colnames(df)
for(n in colnames(df)[-1]){
  p.vals[[n]] <- fisher.test(df[,'de_list'], df[,n],)$p.value
}

go_pvals <- p.vals[grep('^GO_', names(p.vals))]
top_go <- names(head(unlist(go_pvals)[order(unlist(go_pvals))], 100)) # top 100 go terms

np_pvals <- p.vals[grep('^NetPath_', names(p.vals))]
good_np <- names(unlist(np_pvals)[unlist(np_pvals) < 0.01])

# compute overlap / proximity ---------------------------------------------------
inters_mat <- matrix(nrow = length(top_go) + length(good_np), 
                     ncol = length(top_go) + length(good_np), 
                     dimnames = list(c(top_go, good_np), c(top_go, good_np)))

union_mat <- matrix(nrow = length(top_go) + length(good_np), 
                    ncol = length(top_go) + length(good_np), 
                    dimnames = list(c(top_go, good_np), c(top_go, good_np)))

intersection_base <- c(go_terms[top_go], netpath[good_np])

for(row in rownames(inters_mat)){
  for(col in colnames(inters_mat)){
    # print(length(intersect(intersection_base[row], intersection_base[col])))
    inters_mat[row, col] <- length(intersect(intersection_base[[row]], intersection_base[[col]]))
    union_mat[row, col] <- length(union(intersection_base[[row]], intersection_base[[col]]))
  }
}

proximity <- (inters_mat/union_mat)


# network embedding -------------------------------------------------------
# given the low overlap of GO terms and NetPath, in order to embed netpath terms in the network, I need to use a lower threshold for these terms
# higher threshold for GO terms (higher redundancy requires more filtering to display the underlying structure)

# initial threshold of 0.05
adj_mat <- proximity
diag(adj_mat) <- 0
GO_rows <- grep('^GO_', rownames(adj_mat), value = T)
netpath_rows <- grep('^NetPath', rownames(adj_mat), value = T)

# GO vs GO: 0.25
tmp <- adj_mat[GO_rows, GO_rows]
tmp[tmp < 0.25] <- 0
adj_mat[GO_rows, GO_rows] <- tmp

# NetPath vs NetPath: 0.1
tmp <- adj_mat[netpath_rows, netpath_rows]
tmp[tmp < 0.1] <- 0
adj_mat[netpath_rows, netpath_rows] <- tmp

# GO vs NetPath: top 5 links
tmp <- adj_mat[GO_rows, netpath_rows]
for(np in netpath_rows){
  thresh <- sort(adj_mat[GO_rows,np],decreasing = T)[5]
  newcol <- tmp[,np]
  newcol[newcol < thresh] <- 0
  tmp[,np] <- newcol
}

adj_mat[GO_rows, netpath_rows] <- tmp
adj_mat[netpath_rows, GO_rows] <- t(tmp)
adj_mat[adj_mat>0] <- 1

net <- graph_from_adjacency_matrix(adj_mat,mode = 'undirected')

# annotate pathways by database of origin
tmp <- gsub('NetPath.*', 'NetPath',get.vertex.attribute(net, 'name'))
tmp <- gsub('GO_.*','GO', tmp)
net <- set_vertex_attr(net, 'source', value = tmp)

# simplify names
new_names <- get.vertex.attribute(net,'name')
new_names <- gsub('^GO_|^NetPath.','',new_names)
net <- set_vertex_attr(net,'name2',value = new_names)

# visualize on Cytoscape
library(RCy3)
createNetworkFromIgraph(net,"myIgraph")
