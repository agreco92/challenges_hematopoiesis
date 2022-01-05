library(scran)
library(tidyverse)
# differential expression testing 
sce_d1 <- readRDS('./data/day_1_processed/sce2.rds')
sce_d3 <- readRDS('./data/day_3_processed/sce2.rds')


# day 1 -------------------------------------------------------------------
## DEG test
hsc_d1 <- findMarkers(sce_d1[,sce_d1$population == 'HSC'], groups = sce_d1[,sce_d1$population == 'HSC']$condition, test.type = 'wilcox', row.data = rowData(sce_d1))
st_d1 <- findMarkers(sce_d1[,sce_d1$population == 'ST'], groups = sce_d1[,sce_d1$population == 'ST']$condition, test.type = 'wilcox', row.data = rowData(sce_d1))
mpp_d1 <- findMarkers(sce_d1[,sce_d1$population == 'MPP'], groups = sce_d1[,sce_d1$population == 'MPP']$condition, test.type = 'wilcox', row.data = rowData(sce_d1))
cmp_mk_d1 <- findMarkers(sce_d1[,sce_d1$population2 == 'CMP_mkery'], groups = sce_d1[,sce_d1$population2 == 'CMP_mkery']$condition, test.type = 'wilcox', row.data = rowData(sce_d1))
cmp_my_d1 <- findMarkers(sce_d1[,sce_d1$population2 == 'CMP_my'], groups = sce_d1[,sce_d1$population2 == 'CMP_my']$condition ,test.type = 'wilcox', row.data = rowData(sce_d1))

## upregulation
d1_hsc_list <- hsc_d1$sepsis[(hsc_d1$sepsis$FDR <1e-5 & abs(hsc_d1$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d1_st_list <- st_d1$sepsis[(st_d1$sepsis$FDR <1e-5 & abs(st_d1$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d1_mpp_list <- mpp_d1$sepsis[(mpp_d1$sepsis$FDR <1e-5 & abs(mpp_d1$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d1_cmp_mk_list <- cmp_mk_d1$sepsis[(cmp_mk_d1$sepsis$FDR <1e-5 & abs(cmp_mk_d1$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d1_cmp_my_list <- cmp_my_d1$sepsis[(cmp_my_d1$sepsis$FDR <1e-5 & abs(cmp_my_d1$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]

d1_hsc_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/hsc_d1_up.csv')
d1_st_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/st_d1_up.csv')
d1_mpp_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/mpp_d1_up.csv')
d1_cmp_my_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_my_d1_up.csv')
d1_cmp_mk_list %>% as.data.frame()%>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_mk_d1_up.csv')


### list of genes simultaneously upregulated in hsc,st,mpp and myeloid cmps (day 1)
Reduce(intersect, list(d1_hsc_list$Symbol, d1_st_list$Symbol, d1_mpp_list$Symbol, d1_cmp_my_list$Symbol)) %>% 
  as.data.frame() %>% write_csv('./docs/DEGs/hsc_st_mpp_cmp.my_shared_d1_up.csv',col_names = F)

## downregulation
d1_hsc_list <- hsc_d1$sepsis[(hsc_d1$sepsis$FDR <1e-5 & abs(hsc_d1$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d1_st_list <- st_d1$sepsis[(st_d1$sepsis$FDR <1e-5 & abs(st_d1$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d1_mpp_list <- mpp_d1$sepsis[(mpp_d1$sepsis$FDR <1e-5 & abs(mpp_d1$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d1_cmp_mk_list <- cmp_mk_d1$sepsis[(cmp_mk_d1$sepsis$FDR <1e-5 & abs(cmp_mk_d1$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d1_cmp_my_list <- cmp_my_d1$sepsis[(cmp_my_d1$sepsis$FDR <1e-5 & abs(cmp_my_d1$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]

d1_hsc_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/hsc_d1_dn.csv')
d1_st_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/st_d1_dn.csv')
d1_mpp_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/mpp_d1_dn.csv')
d1_cmp_my_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_my_d1_dn.csv')
d1_cmp_mk_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_mk_d1_dn.csv')


# day 3 -------------------------------------------------------------------

## DEG test
hsc_d3 <- findMarkers(sce_d3[,sce_d3$population == 'HSC'], groups = sce_d3[,sce_d3$population == 'HSC']$condition, test.type = 'wilcox', row.data = rowData(sce_d3))
st_d3 <- findMarkers(sce_d3[,sce_d3$population == 'ST'], groups = sce_d3[,sce_d3$population == 'ST']$condition, test.type = 'wilcox', row.data = rowData(sce_d3))
mpp_d3 <- findMarkers(sce_d3[,sce_d3$population == 'MPP'], groups = sce_d3[,sce_d3$population == 'MPP']$condition, test.type = 'wilcox', row.data = rowData(sce_d3))
cmp_mk_d3 <- findMarkers(sce_d3[,sce_d3$population2 == 'CMP_mkery'], groups = sce_d3[,sce_d3$population2 == 'CMP_mkery']$condition, test.type = 'wilcox', row.data = rowData(sce_d3))
cmp_my_d3 <- findMarkers(sce_d3[,sce_d3$population2 == 'CMP_my'], groups = sce_d3[,sce_d3$population2 == 'CMP_my']$condition ,test.type = 'wilcox', row.data = rowData(sce_d3))

## upregulation
d3_hsc_list <- hsc_d3$sepsis[(hsc_d3$sepsis$FDR <1e-5 & abs(hsc_d3$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d3_st_list <- st_d3$sepsis[(st_d3$sepsis$FDR <1e-5 & abs(st_d3$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d3_mpp_list <- mpp_d3$sepsis[(mpp_d3$sepsis$FDR <1e-5 & abs(mpp_d3$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d3_cmp_mk_list <- cmp_mk_d3$sepsis[(cmp_mk_d3$sepsis$FDR <1e-5 & abs(cmp_mk_d3$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]
d3_cmp_my_list <- cmp_my_d3$sepsis[(cmp_my_d3$sepsis$FDR <1e-5 & abs(cmp_my_d3$sepsis$AUC.ctrl) > 2/3), c('Symbol','FDR','AUC.ctrl')]

d3_hsc_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/hsc_d3_up.csv')
d3_st_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/st_d3_up.csv')
d3_mpp_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/mpp_d3_up.csv')
d3_cmp_my_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_my_d3_up.csv')
d3_cmp_mk_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_mk_d3_up.csv')

## downregulation
d3_hsc_list <- hsc_d3$sepsis[(hsc_d3$sepsis$FDR <1e-5 & abs(hsc_d3$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d3_st_list <- st_d3$sepsis[(st_d3$sepsis$FDR <1e-5 & abs(st_d3$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d3_mpp_list <- mpp_d3$sepsis[(mpp_d3$sepsis$FDR <1e-5 & abs(mpp_d3$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d3_cmp_mk_list <- cmp_mk_d3$sepsis[(cmp_mk_d3$sepsis$FDR <1e-5 & abs(cmp_mk_d3$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]
d3_cmp_my_list <- cmp_my_d3$sepsis[(cmp_my_d3$sepsis$FDR <1e-5 & abs(cmp_my_d3$sepsis$AUC.ctrl) < 1/3), c('Symbol','FDR','AUC.ctrl')]

d3_hsc_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/hsc_d3_dn.csv')
d3_st_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/st_d3_dn.csv')
d3_mpp_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/mpp_d3_dn.csv')
d3_cmp_my_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_my_d3_dn.csv')
d3_cmp_mk_list %>% as.data.frame() %>% format(digits = 3) %>% write_csv('./docs/DEGs/cmp_mk_d3_dn.csv')