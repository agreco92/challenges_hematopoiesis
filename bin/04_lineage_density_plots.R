library(Seurat)
library(scater)
library(scran)
library(patchwork)
lineage.modules = fgsea::gmtPathways(gmt.file = './data/lin_sets/lineage.sets.gmt')

# day ^1 ------------------------------------------------------------------
ctrl <- readRDS('./data/day_1_processed/ctrl.rds')
sepsis <- readRDS('./data/day_1_processed/sepsis.rds')
sce_d1 <- readRDS('./data/day_1_processed/sce2.rds')

# compute lineage scores
seu <- as.Seurat(sce_d1)
for(i in names(lineage.modules)){
  seu <- AddModuleScore(seu, features = list(lineage.modules[[i]]), name = i)
}
colData(sce_d1)[,names(lineage.modules)] <- seu@meta.data[,paste0(names(lineage.modules),'1')]

## scale score by their variance
#d1_scale <- apply(colData(sce_d1)[sce_d1$condition == 'ctrl',names(lineage.modules)], MARGIN = 2, FUN = sd)
d1_scale <- sd(unlist(colData(sce_d1)[sce_d1$condition == 'ctrl',names(lineage.modules)]))

tmp <- t(t(as.matrix(colData(sce_d1)[,names(lineage.modules)]))/d1_scale)
colData(sce_d1)[,names(lineage.modules)] <- tmp

df <- data.frame(condition = sce_d1$condition, population = sce_d1$population2) %>% 
  bind_cols(as.data.frame(colData(sce_d1)[,names(lineage.modules)]))
df2 <- df %>% gather(stem, ly, my, ery, mk, mep, key = 'module', value = 'score')
df3 <- df2 %>% group_by(condition, population, module) %>% summarize(score = median(score)) %>% ungroup()
df4 <- df3 %>% spread(condition, score) %>% mutate(sepsis = sepsis - ctrl, ctrl = 0) %>%  gather(sepsis, ctrl, key = 'condition', value = 'score')
df4$population <- as.character(df4$population)

# density plot (day1)
semi_void_theme <- theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank())

# myeloid row
hsc_my <- df %>% dplyr::filter(population == 'HSC') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
st_my <- df %>% dplyr::filter(population == 'ST') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
mpp_my <- df %>% dplyr::filter(population == 'MPP') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.my_my <- df %>% dplyr::filter(population == 'CMP_my') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.mkery_my <- df %>% dplyr::filter(population == 'CMP_mkery') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme #+ guides(color =F)

myeloid_row <- (hsc_my + st_my + mpp_my + cmp.my_my + cmp.mkery_my) + plot_layout(ncol =5) 

# mkery row
hsc_mkery <- df %>% dplyr::filter(population == 'HSC') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
st_mkery <- df %>% dplyr::filter(population == 'ST') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
mpp_mkery <- df %>% dplyr::filter(population == 'MPP') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.my_mkery <- df %>% dplyr::filter(population == 'CMP_my') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.mkery_mkery <- df %>% dplyr::filter(population == 'CMP_mkery') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme #+ guides(color =F)

mkery_row <- (hsc_mkery + st_mkery + mpp_mkery + cmp.my_mkery + cmp.mkery_mkery) + plot_layout(ncol =5) 

# lymphoid row
hsc_ly <- df %>% dplyr::filter(population == 'HSC') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
st_ly <- df %>% dplyr::filter(population == 'ST') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
mpp_ly <- df %>% dplyr::filter(population == 'MPP') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.my_ly <- df %>% dplyr::filter(population == 'CMP_my') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.mkery_ly <- df %>% dplyr::filter(population == 'CMP_mkery') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme #+ guides(color =F)

ly_row <- (hsc_ly + st_ly + mpp_ly + cmp.my_ly + cmp.mkery_ly) + plot_layout(ncol =5) 

# compose output
density_plot_d1 <- (myeloid_row / mkery_row / ly_row) 

# day 3 -------------------------------------------------------------------
ctrl <- readRDS('./data/day_3_processed/ctrl.rds')
sepsis <- readRDS('./data/day_3_processed/sepsis.rds')
sce_d3 <- readRDS('./data/day_3_processed/sce2.rds')

# compute lineage scores
seu <- as.Seurat(sce_d3)
for(i in names(lineage.modules)){
  seu <- AddModuleScore(seu, features = list(lineage.modules[[i]]), name = i)
}
colData(sce_d3)[,names(lineage.modules)] <- seu@meta.data[,paste0(names(lineage.modules),'1')]

## scale score by their variance
d3_scale <- sd(unlist(colData(sce_d3)[sce_d1$condition == 'ctrl',names(lineage.modules)]))
tmp <- t(t(as.matrix(colData(sce_d3)[,names(lineage.modules)]))/d3_scale)
colData(sce_d3)[,names(lineage.modules)] <- tmp

df <- data.frame(condition = sce_d3$condition, population = sce_d3$population2) %>% 
  bind_cols(as.data.frame(colData(sce_d3)[,names(lineage.modules)]))
df2 <- df %>% gather(stem, ly, my, ery, mk, mep, key = 'module', value = 'score')
df3 <- df2 %>% group_by(condition, population, module) %>% summarize(score = median(score)) %>% ungroup()
df4 <- df3 %>% spread(condition, score) %>% mutate(sepsis = sepsis - ctrl, ctrl = 0) %>%  gather(sepsis, ctrl, key = 'condition', value = 'score')
df4$population <- as.character(df4$population)

# density plot (day1)
semi_void_theme <- theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank())

# myeloid row
hsc_my <- df %>% dplyr::filter(population == 'HSC') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
st_my <- df %>% dplyr::filter(population == 'ST') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
mpp_my <- df %>% dplyr::filter(population == 'MPP') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.my_my <- df %>% dplyr::filter(population == 'CMP_my') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.mkery_my <- df %>% dplyr::filter(population == 'CMP_mkery') %>% ggplot(aes(x = my, col = condition)) + geom_density() + semi_void_theme #+ guides(color =F)

myeloid_row <- (hsc_my + st_my + mpp_my + cmp.my_my + cmp.mkery_my) + plot_layout(ncol =5) 

# mkery row
hsc_mkery <- df %>% dplyr::filter(population == 'HSC') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
st_mkery <- df %>% dplyr::filter(population == 'ST') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
mpp_mkery <- df %>% dplyr::filter(population == 'MPP') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.my_mkery <- df %>% dplyr::filter(population == 'CMP_my') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.mkery_mkery <- df %>% dplyr::filter(population == 'CMP_mkery') %>% ggplot(aes(x = mep, col = condition)) + geom_density() + semi_void_theme #+ guides(color =F)

mkery_row <- (hsc_mkery + st_mkery + mpp_mkery + cmp.my_mkery + cmp.mkery_mkery) + plot_layout(ncol =5) 

# lymphoid row
hsc_ly <- df %>% dplyr::filter(population == 'HSC') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
st_ly <- df %>% dplyr::filter(population == 'ST') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
mpp_ly <- df %>% dplyr::filter(population == 'MPP') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.my_ly <- df %>% dplyr::filter(population == 'CMP_my') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme + guides(color =F)
cmp.mkery_ly <- df %>% dplyr::filter(population == 'CMP_mkery') %>% ggplot(aes(x = ly, col = condition)) + geom_density() + semi_void_theme #+ guides(color =F)

ly_row <- (hsc_ly + st_ly + mpp_ly + cmp.my_ly + cmp.mkery_ly) + plot_layout(ncol =5) 

# compose output
density_plot_d3 <- (myeloid_row / mkery_row / ly_row) 



