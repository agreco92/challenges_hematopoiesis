library(scater)
ctrl <- readRDS('./data/day_3_processed/ctrl.rds')
sepsis <- readRDS('./data/day_3_processed/sepsis.rds')
sce <- readRDS('./data/day_3_processed/sce.rds') # cells in both conditions

# merge cmp classification from ctrl and sepsis dataset
## assign -1 (unknown branch) to all cells
sce$mk_gmp_branch  <- '-1'

## assign consistent cell names
colnames(ctrl) <- paste0(ctrl$Sample,'_', ctrl$Barcode)
colnames(sepsis) <- paste0(sepsis$Sample,'_', sepsis$Barcode)
colnames(sce) <- paste0(sce$Sample, '_',sce$Barcode)
sce <- sce[,c(colnames(ctrl),colnames(sepsis))]

## merge ctrl and sepsis branch classification into merged object
colData(sce)[,'mk_gmp_branch'] <- c(ctrl$mk_gmp_branch,sepsis$mk_gmp_branch)

# subdivide CMP based on branch
sce$population2 <- 'CMP'
sce$population2[sce$population == 'HSC'] = 'HSC'
sce$population2[sce$population == 'ST'] = 'ST'
sce$population2[sce$population == 'MPP'] = 'MPP'

sce$population2[(sce$population == 'CMP' & sce$mk_gmp_branch == '2') & sce$condition == 'ctrl' ] = 'CMP_mkery'
sce$population2[(sce$population == 'CMP' & sce$mk_gmp_branch == '2') & sce$condition == 'sepsis' ] = 'CMP_mkery'

sce$population2[(sce$population == 'CMP' & sce$mk_gmp_branch == '3') & sce$condition == 'ctrl'] = 'CMP_my'
sce$population2[(sce$population == 'CMP' & sce$mk_gmp_branch == '3') & sce$condition == 'sepsis'] = 'CMP_my'

plotDiffusionMap(sce, col = 'population2')

# save Single Cell Experiment object
saveRDS(sce, './data/day_3_processed/sce2.rds')
