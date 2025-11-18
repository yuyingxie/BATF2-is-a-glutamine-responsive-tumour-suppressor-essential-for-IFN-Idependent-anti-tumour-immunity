library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(matrixStats)
library(parallel)
library(reshape2)
library(ComplexHeatmap)
library(slingshot)
library(qlcMatrix)
set.seed(1)

library(phateR) 
library(dplyr)   # define funciton %>%
library(hdf5r)
 X11(type='cairo')

dat_saver = readRDS('batf.denoised_SAVERX.rds')
dat = dat_saver$estimate
gene = readRDS('gene_name.rds') # this is the gene list from Batf2_merge.rds
dat = dat[gene, ]
saveRDS(dat, file = 'Batf2_saver.rds')

dat_saver = readRDS('Batf2_saver.rds')
MOC1_cells = readRDS('MOC1_cell_name.rds')
MOC2_cells = readRDS('MOC2_cell_name.rds')
KO1_cells = readRDS('KO1_cell_name.rds')
KO2_cells = readRDS('KO2_cell_name.rds')

id = match(substr(MOC1_cells, start = 6, stop = 21), substr(colnames(dat_saver), 1, 16))
colnames(dat_saver)[id]  # -1 is MOC1
id = match(substr(MOC2_cells, start = 6, stop = 21), substr(colnames(dat_saver), 1, 16))
colnames(dat_saver)[id]  # -2 is MOC2
id = match(substr(KO1_cells, start = 7, stop = 22), substr(colnames(dat_saver), 1, 16))
colnames(dat_saver)[id]  # -3 is KO1
id = match(substr(KO2_cells, start = 7, stop = 22), substr(colnames(dat_saver), 1, 16))
colnames(dat_saver)[id]  # -3 is KO1

id = substr(colnames(dat_saver), 17, 18) == '-1'
MOC1 = dat_saver[, id]

dat_merge = readRDS('Batf2.merge.rds')
dat.list <- SplitObject(dat_merge, split.by = "orig.ident")
A1 = dat.list[[1]]
dat1 = GetAssayData(A1, slot = 'counts')
id = match(substr(colnames(dat1), 6, 21), substr(colnames(MOC1), 1, 16))
plot(dat1[1, ], MOC1[1, id])


dat <- CreateSeuratObject(counts = dat, min.cells = 10)
dat[['percent.mt']] <- PercentageFeatureSet(dat, pattern = "^mt-")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Sample_id = rep('MOC1', dim(dat)[2])
Treatment = rep('MOC', dim(dat)[2])
meta = data.frame(Sample = Sample_id,  Treatment = Treatment,
                  row.names = colnames(dat), Experiment = rep('Batf2', dim(dat)[2] ))
dat = AddMetaData(dat, metadata = meta)
dat.list[[sample_name]] = dat


##############################
## From Count
#########################

dat_merge = readRDS('Batf2.merge.rds')
dat.list <- SplitObject(dat_merge, split.by = "orig.ident")
A1 = dat.list[[1]]
saveRDS(colnames(A1), file = 'MOC1_cell_name.rds')
A2 = dat.list[[2]]
saveRDS(colnames(A2), file = 'MOC2_cell_name.rds')
A3 = dat.list[[3]]
saveRDS(colnames(A3), file = 'KO1_cell_name.rds')
A4 = dat.list[[4]]
saveRDS(colnames(A4), file = 'KO2_cell_name.rds')



dat_list = list()
dat_list[['MOC1']] = dat.list[[1]]
dat_list[['MOC2']] = dat.list[[2]]
dat_list[['Batf2_1']] = dat.list[[3]]
dat_list[['Batf2_2']] = dat.list[[4]]
rm(dat.list)
DefaultAssay(dat_list[[1]] )


dat.list = list()
i = 1
sample.tmp.seurat = dat_list[[i]]
sample_name = names(dat_list)[i]
print(sample_name)
dat = GetAssayData(sample.tmp.seurat, slot = 'counts')
dat <- CreateSeuratObject(counts = dat, min.cells = 10)
dat[['percent.mt']] <- PercentageFeatureSet(dat, pattern = "^mt-")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Sample_id = rep('MOC1', dim(dat)[2])
Treatment = rep('MOC', dim(dat)[2])
meta = data.frame(Sample = Sample_id,  Treatment = Treatment,
                  row.names = colnames(dat), Experiment = rep('Batf2', dim(dat)[2] ))
dat = AddMetaData(dat, metadata = meta)
dat.list[[sample_name]] = dat

i = 2
sample.tmp.seurat = dat_list[[i]]
sample_name = names(dat_list)[i]
print(sample_name)
dat = GetAssayData(sample.tmp.seurat, slot = 'counts')
dat <- CreateSeuratObject(counts = dat, min.cells = 10)
dat[['percent.mt']] <- PercentageFeatureSet(dat, pattern = "^mt-")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Sample_id = rep('MOC2', dim(dat)[2])
Treatment = rep('MOC', dim(dat)[2])
meta = data.frame(Sample = Sample_id,  Treatment = Treatment,
                  row.names = colnames(dat), Experiment = rep('Batf2', dim(dat)[2] ))
dat = AddMetaData(dat, metadata = meta)
dat.list[[sample_name]] = dat

i = 3
sample.tmp.seurat = dat_list[[i]]
sample_name = names(dat_list)[i]
print(sample_name)
dat = GetAssayData(sample.tmp.seurat, slot = 'counts')
dat <- CreateSeuratObject(counts = dat, min.cells = 10)
dat[['percent.mt']] <- PercentageFeatureSet(dat, pattern = "^mt-")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(dat$nCount_RNA, dat$nFeature_RNA)
abline(v = 65000, col = 'red')
id = dat$nCount_RNA > 65000 
dat = subset(x = dat, cells = colnames(dat)[!id]  )
FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Sample_id = rep('KO1', dim(dat)[2])
Treatment = rep('KO', dim(dat)[2])
meta = data.frame(Sample = Sample_id,  Treatment = Treatment,
                  row.names = colnames(dat), Experiment = rep('Batf2', dim(dat)[2] ))
dat = AddMetaData(dat, metadata = meta)
dat.list[[sample_name]] = dat

i = 4
sample.tmp.seurat = dat_list[[i]]
sample_name = names(dat_list)[i]
print(sample_name)
dat = GetAssayData(sample.tmp.seurat, slot = 'counts')
dat <- CreateSeuratObject(counts = dat, min.cells = 10)
dat[['percent.mt']] <- PercentageFeatureSet(dat, pattern = "^mt-")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Sample_id = rep('KO2', dim(dat)[2])
Treatment = rep('KO', dim(dat)[2])
meta = data.frame(Sample = Sample_id,  Treatment = Treatment,
                  row.names = colnames(dat), Experiment = rep('Batf2', dim(dat)[2] ))
dat = AddMetaData(dat, metadata = meta)
dat.list[[sample_name]] = dat
rm(dat_list)


############################################
###### Cell cylce
#########################################

for(i in 1:length(dat.list)){
    sample.tmp.seurat = dat.list[[i]]
    sample_name = names(dat.list)[i]
    print(sample_name)
    sample.tmp.seurat <- SCTransform(sample.tmp.seurat, vars.to.regress = "percent.mt", verbose = TRUE,
     variable.features.n = 5000)
    dat.list[sample_name] = sample.tmp.seurat
}

saveRDS(dat.list, file = 'dat_list_SCT_03_10_2021.rds')
dat.list = readRDS('dat_list_SCT_03_10_2021.rds')

############################################################
### Integration  ###########################################
############################################################

dat.features <- SelectIntegrationFeatures(object.list = dat.list, nfeatures = 3000)

gene_list = intersect(mgene, rownames(dat.list[[1]])) 
for(i in 2:4){
    gene_list = intersect(gene_list, rownames(dat.list[[i]]))    
}

dat.features = unique(c(gene_list, dat.features))

dat.list <- PrepSCTIntegration(object.list = dat.list, anchor.features = dat.features)

#Anchors and Integration
dat.anchors <- FindIntegrationAnchors(object.list = dat.list, normalization.method = "SCT", 
                                        anchor.features = dat.features)
dat.integrated <- IntegrateData(anchorset = dat.anchors, normalization.method = "SCT")

rm(dat.anchors)
rm(dat.list)

dat.integrated <- RunPCA(object = dat.integrated, verbose = FALSE, npcs = 100)
dat.integrated <- RunUMAP(object = dat.integrated, dims = 1:50)

###cluster
dat.integrated <- FindNeighbors(object = dat.integrated, dims = 1:50)
dat.integrated <- FindClusters(object = dat.integrated, resolution = 0.8)### 

saveRDS(dat.integrated, file = 'dat_integrated_03_09_2021.rds')
dat.integrated = readRDS('dat_integrated_03_09_2021.rds')



pdf('Integrated_sample_03_09_2021.pdf')
DimPlot(dat.integrated, group.by = c("Sample"), combine = FALSE)
dev.off()

pdf('Integrated_Treatment_03_20_2021.pdf')
DimPlot(dat.integrated, group.by =  "Treatment")
dev.off()

FeaturePlot(dat.integrated, features = c("Ifnar1"))


DimPlot(dat.integrated, group.by = c("state"), combine = FALSE)

DefaultAssay(dat.integrated) = 'integrated'

pdf('Integrated_cluster_03_09_2021.pdf')
DimPlot(dat.integrated, label = TRUE)
dev.off()

#########################################################################
##########  Annotation 
##################################################################
library(SingleR)
DefaultAssay(dat.integrated) = 'RNA'
ImmGe = ImmGenData() # This is for mouse 
ref <- ImmGe
logdata <- GetAssayData(dat.integrated[['SCT']], slot = 'data')
dat_clusters = dat.integrated$seurat_clusters

singler_cluster <- SingleR(test = logdata, ref = ref, labels = ref$label.main,  
                           method = "cluster", clusters = dat_clusters, de.method="wilcox")

save(singler_cluster,file = 'SingleR_Label_pred.rda')
write.csv(singler_cluster,'singler_pred_labels_03_11_2021.csv')

common <- intersect(rownames(ref), rownames(dat.integrated))
logdata_sub <- logdata[common,]
ref <- ref[common, ]
singler_cluster2 <- SingleR(test = logdata_sub, ref = ref, labels = ref$label.main,  
                           method = "cluster", clusters = T_clusters, de.method="wilcox")
save(singler_cluster2,file = 'SingleR_2_Label_pred.rda')

write.csv(singler_cluster2, file = 'singler_pred_label2_11_27_2020.csv')
###################################
## Find Marker
##################################
DefaultAssay(dat.integrated) <- "RNA"
print('***tsne and umap***')

dat.integrated@misc$markers <- FindAllMarkers(object = dat.integrated, assay = 'RNA', only.pos = TRUE, test.use = 'wilcox')
write.table(dat.integrated@misc$markers,file='marker_wilcox.txt',row.names = FALSE,quote = FALSE,sep = '\t')

dpi = 300
png(file="feature_vln_03_11_2020.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = dat.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

hc.markers = read.delim2("marker_wilcox.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50
tt1 = DoHeatmap(object = subset(dat.integrated, downsample = 500), features = top50$gene) + NoLegend()
ggplot2::ggsave(file="marker_heatmap_MAST.pdf", plot = tt1, device = 'pdf',width = 20,
           height = 16, units = "in",dpi = dpi,limitsize = FALSE)
print('***find markers***')

top_gene = hc.markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_log2FC)
write.csv(cbind(top_gene$gene, top_gene$cluster), file = 'genes_list_Ifnar1_03_11_2021.csv')

res =Get_mean(dat.integrated)
write.csv(res, 'Ifnar1_mean_sd_03_11_2021.csv')


##############################################################
#######   Subset Cluster 14
cluster_id = 13
id = dat.integrated$seurat_clusters %in% cluster_id
dat_sub = subset(dat.integrated, cells = colnames(dat.integrated)[id])

VlnPlot(dat_sub, features = 'Ifnar1', assay = 'RNA', group.by = 'Treatment',  pt.size = 0.02)

DefaultAssay(dat_sub) = 'RNA'

VlnPlot(dat_sub, features = 'Cxcl9',  group.by = 'Sample',  pt.size = 0.02)

DefaultAssay(dat_sub) = 'RNA'
VlnPlot(dat_sub, features = 'Isg15',  group.by = 'Sample',  pt.size = 0.02)

FeaturePlot(dat.integrated, features = c("Ifnar1", "Cxcl9", "Cxcl10", "Oasl1"))



FeaturePlot(dat.integrated, features = c(""))




###########################################################
##### Composition Analysis
###########################################################

library(compositions)

B = table(dat.integrated$Treatment, dat.integrated$seurat_clusters)
B1 = table(dat.integrated$Sample, dat.integrated$seurat_clusters)

B1_sum = apply(B1, 1, sum)
B1_relative = B1 / B1_sum

A = matrix(0, 20, 2)
colnames(A) = c('KO', 'MOC')

for(i in 1:20){
    A[i, 1] = mean(B1_relative[1:2, i])
    A[i, 2] = mean(B1_relative[3:4, i])

}
rownames(A) = 1:20

result = matrix(0, 20, 2)
colnames(result) = c('cluster', 'Pval')

result[, 1] = 1:20

for(i in 1:20){
     result[i, 2] = t.test(B1_relative[1:2, i], B1_relative[3:4, i])$p.value
}


#############################################
##### Bar plot
library(RColorBrewer)

bar.col = c(brewer.pal(12, 'Set3'))


pdf("Cluster_bar_mean_horiz_03_20_2021.pdf", width = 20, height = 7)
opar <- par(lwd = 0.0003)
par(mar = c(5, 6, 4.1, 15), xpd = TRUE)
barplot(height = A, col = bar.col, width = 2, axes = F, xlab = "",  space = 0.1,
      border = NA, las = 2, horiz = T,  names.arg=c("KO", "MOC")) 
#axis(2)
axis(1)
legend(x = 1.05, y = 4, 
		legend = (rownames(A)), #in order from top to bottom
		fill = bar.col, # 6:1 reorders so legend order matches graph
		title = "Cluster", ncol = 2, cex = 1.2, bty = "n")
dev.off()

