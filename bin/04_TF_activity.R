# includes ----------------------------------------------------------------
library(scran)
library(viper)

# import data -------------------------------------------------------------
# Dorothea mouse regulons
regulons <- readRDS('./data/pathways/dorothea_mouse_AB_viper_format.rds')
sce <- readRDS('./data/day_1_processed/sce2.rds')

exp.mat <- scale(logcounts(sce))
rownames(exp.mat) <- rowData(sce)[,'Symbol']

# compute regulons activity using viper
TF_viper <-  viper(eset = exp.mat, regulon = regulons, minsize = 5)
colnames <- paste0('v_',rownames(TF_viper))
colData(sce)[colnames] = t(TF_viper)

df <- colData(sce)[,c(colnames)]

# HSC visualization---------------------------------------------------------------------
ctrl_mask <- sce$condition == 'ctrl' & sce$population == 'HSC'
d1_mask <- sce$condition == 'sepsis' & sce$population == 'HSC'

# differential activity
HSC_pvals <- apply(df, 2, function(x){wilcox.test(x[d1_mask], x[ctrl_mask])$p.value})
HSC_diff <- apply(df, 2, function(x){t.test(x[d1_mask], x[ctrl_mask])$statistic})

names(HSC_pvals) = names(HSC_diff) <- gsub('v_','', names(HSC_pvals))

HSC_pvals <- HSC_pvals %>% as.data.frame() %>% rownames_to_column('tf') %>% mutate(t = HSC_diff)
colnames(HSC_pvals) <- c('tf','p_val','t_stat')

HSC_pvals <- HSC_pvals %>% dplyr::arrange(p_val)
ggplot(HSC_pvals %>% dplyr::filter(p_val < 1e-20), aes(x = reorder(tf, t_stat, sum), y = t_stat)) + geom_col(aes(fill = -log10(1e-300+p_val))) + coord_flip()

#### ST #####
ctrl_mask <- sce$condition == 'ctrl' & sce$population == 'ST'
d1_mask <- sce$condition == 'sepsis' & sce$population == 'ST'

ST_pvals <- apply(df, 2, function(x){wilcox.test(x[d1_mask], x[ctrl_mask])$p.value})
ST_diff <- apply(df, 2, function(x){t.test(x[d1_mask], x[ctrl_mask])$statistic})

names(ST_pvals) = names(ST_diff) <- gsub('v_','', names(ST_pvals))

ST_pvals <- ST_pvals %>% as.data.frame() %>% rownames_to_column('tf') %>% mutate(t = ST_diff)
colnames(ST_pvals) <- c('tf','p_val','t_stat')

ST_pvals <- ST_pvals %>% dplyr::arrange(p_val)
ggplot(ST_pvals %>% dplyr::filter(p_val < 1e-20), aes(x = reorder(tf, t_stat, sum), y = t_stat)) + geom_col(aes(fill = -log10(1e-300+p_val))) + coord_flip()

##### MPP #####
ctrl_mask <- sce$condition == 'ctrl' & sce$population == 'MPP'
d1_mask <- sce$condition == 'sepsis' & sce$population == 'MPP'

MPP_pvals <- apply(df, 2, function(x){wilcox.test(x[d1_mask], x[ctrl_mask])$p.value})
MPP_diff <- apply(df, 2, function(x){t.test(x[d1_mask], x[ctrl_mask])$statistic})

names(MPP_pvals) = names(MPP_diff) <- gsub('v_','', names(MPP_pvals))

MPP_pvals <- MPP_pvals %>% as.data.frame() %>% rownames_to_column('tf') %>% mutate(t = MPP_diff)
colnames(MPP_pvals) <- c('tf','p_val','t_stat')

MPP_pvals <- MPP_pvals %>% dplyr::arrange(p_val)
ggplot(MPP_pvals %>% dplyr::filter(p_val < 1e-20), aes(x = reorder(tf, t_stat, sum), y = t_stat)) + geom_col(aes(fill = -log10(1e-300+p_val))) + coord_flip()

##### CMP_my #####
ctrl_mask <- sce$condition == 'ctrl' & sce$population2 == 'CMP_my'
d1_mask <- sce$condition == 'sepsis' & sce$population2 == 'CMP_my'

CMP_my_pvals <- apply(df, 2, function(x){wilcox.test(x[d1_mask], x[ctrl_mask])$p.value})
CMP_my_diff <- apply(df, 2, function(x){t.test(x[d1_mask], x[ctrl_mask])$statistic})

names(CMP_my_pvals) = names(CMP_my_diff) <- gsub('v_','', names(CMP_my_pvals))

CMP_my_pvals <- CMP_my_pvals %>% as.data.frame() %>% rownames_to_column('tf') %>% mutate(t = CMP_my_diff)
colnames(CMP_my_pvals) <- c('tf','p_val','t_stat')

CMP_my_pvals <- CMP_my_pvals %>% dplyr::arrange(p_val)
ggplot(CMP_my_pvals %>% dplyr::filter(p_val < 1e-20), aes(x = reorder(tf, t_stat, sum), y = t_stat)) + geom_col(aes(fill = -log10(1e-300+p_val))) + coord_flip()

##### CMP_mkery #####
ctrl_mask <- sce$condition == 'ctrl' & sce$population2 == 'CMP_mkery'
d1_mask <- sce$condition == 'sepsis' & sce$population2 == 'CMP_mkery'

CMP_mkery_pvals <- apply(df, 2, function(x){wilcox.test(x[d1_mask], x[ctrl_mask])$p.value})
CMP_mkery_diff <- apply(df, 2, function(x){t.test(x[d1_mask], x[ctrl_mask])$statistic})
names(CMP_mkery_pvals) = names(CMP_mkery_diff) <- gsub('v_','', names(CMP_mkery_pvals))

CMP_mkery_pvals <- CMP_mkery_pvals %>% as.data.frame() %>% rownames_to_column('tf') %>% mutate(t = CMP_mkery_diff)
colnames(CMP_mkery_pvals) <- c('tf','p_val','t_stat')

CMP_mkery_pvals <- CMP_mkery_pvals %>% dplyr::arrange(p_val)
ggplot(CMP_mkery_pvals %>% dplyr::filter(p_val < 1e-20), aes(x = reorder(tf, t_stat, sum), y = t_stat)) + geom_col(aes(fill = -log10(1e-300+p_val))) + coord_flip()

# visualize top terms in heatmap (paper figure)-------------------------------------------
## merge analysis on single populations
dff <- rbind(HSC_pvals %>% mutate(population = 'HSC') %>% mutate(rank_p_val = rank(p_val,ties.method = 'min')),
             ST_pvals %>% mutate(population = 'ST') %>% mutate(rank_p_val = rank(p_val,ties.method = 'min')),
             MPP_pvals %>% mutate(population = 'MPP') %>% mutate(rank_p_val = rank(p_val,ties.method = 'min')),
             CMP_my_pvals %>% mutate(population = 'CMP_my') %>% mutate(rank_p_val = rank(p_val,ties.method = 'min')),
             CMP_mkery_pvals %>% mutate(population = 'CMP_mkery') %>% mutate(rank_p_val = rank(p_val,ties.method = 'min'))) %>% as.data.frame()
dff$population = factor(dff$population, levels = c('HSC','ST','MPP','CMP_my','CMP_mkery'))

## select top significant terms
tf_sel <- unique(dff$tf[dff$rank_p_val<=15])

dfff <- dff %>% dplyr::filter(tf%in% tf_sel) %>% dplyr::select(tf, t_stat, population) %>% 
  spread(tf, t_stat) %>% column_to_rownames('population')

dfff <- sign(dfff) * log1p(abs(dfff)) # transform difference in logscale

pheatmap::pheatmap(t(dfff), cluster_cols = F, border_color = NA, fontsize = 8, legend = F)
