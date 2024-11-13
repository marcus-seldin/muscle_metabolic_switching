setwd('G:/My Drive/lab files/Lauren Albrecht/muslce pfk story/git analyses/')
library(SeuratDisk)
library(Seurat)
library(WGCNA)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)
Convert("./raw_data/TS_Muscle.h5ad", ".h5seurat")
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat(".h5seurat", assays = "RNA")

gene1 = 'PFKM'
tissue_name = "Muscle"
pbmc = seuratObject
pdf(file = paste0('all cells UMAP GTEx ',tissue_name, '.pdf'))
DimPlot(pbmc, reduction = "umap",  group.by = 'free_annotation',label = TRUE,
        repel = TRUE)
dev.off()

pdf(file = paste0(gene1, ' expression UMAP GTEx ', tissue_name, '.pdf'))
FeaturePlot(pbmc, features = c(gene1),  reduction = "umap",)
dev.off()
Idents(pbmc) = pbmc@meta.data$free_annotation
table(Idents(pbmc))
pdf(file = paste0(gene1, ' expression Violin GTEx', tissue_name, '.pdf'))
VlnPlot(object = pbmc, features = c(gene1), group.by = 'free_annotation') + theme(legend.position = 'none')
dev.off()

#heatmap for cell cors
sc_matrix = as.data.frame(pbmc[["RNA"]]@counts)
colnames(sc_matrix) = pbmc@meta.data$free_annotation
sc_matrix= as.data.frame(sc_matrix)


#make pathway heatmap
df1 = as.data.frame(c('G6PDH', 'H6PD', 'PGD', 'RPE', 'RPI', 'TKT', 'TALD', 'PGLS', 'PGM1', 'PGM2', 'PRPS1'))
colnames(df1) = 'gene_symbol'
df1$pathways = paste0('PPP')

df2 = as.data.frame(c('GAPDH', 'PHGDH', 'PKM', 'ENO', 'PGK', 'PGM', 'PK'))
colnames(df2) = 'gene_symbol'
df2$pathways = paste0('Glycolysis')

df3 = as.data.frame(c('CS', 'ACO', 'IDH', 'MDH2', 'FH', 'OGDH', 'PDHA', 'SDHC', 'SUCLG'))
colnames(df3) = 'gene_symbol'
df3$pathways = paste0('TCA')

full_paths1 = as.data.frame(rbind(df1, df2, df3))
cor_genes= c('PFKM', 'LAMP1')
cc3 = sc_matrix[row.names(sc_matrix) %in% cor_genes,]

gene_heat1 = sc_matrix[row.names(sc_matrix) %in% full_paths1$gene_symbol,]
cor_incell = function(celt12){
  cc1 = gene_heat1[,colnames(gene_heat) %in% celt12]
  cc2 = as.data.frame(t(cc1))
  
  cc4 = cc3[,colnames(cc3) %in% celt12]
  cc5 = as.data.frame(t(cc4))
  cc6 = bicorAndPvalue(cc2, cc5, use = 'p')
  nn1 = melt(cc6$bicor)
  colnames(nn1) = c('gene_1', 'gene_2', 'bicor')
  nn1$pvalue = melt(cc6$p)$value
  nn1$cellType = paste0(celt12)
  return(nn1)
}

full_heat_tablemelt = as.data.frame(rbind(cor_incell('skeletal muscle satellite stem cell'), cor_incell('slow muscle cell'), cor_incell('fast muscle cell')))
head(full_heat_tablemelt)
full_heat_tablemelt$pathway = full_paths1$pathways[match(full_heat_tablemelt$gene_1, full_paths1$gene_symbol)]

full_heat_tablemelt$gP1 = paste0(full_heat_tablemelt$gene_2, '_', full_heat_tablemelt$cellType)
heat1 = dcast(full_heat_tablemelt,  gene_1 ~ gP1, fun.aggregate = mean, value.var = 'bicor')
row.names(heat1) = heat1$gene_1
heat1$gene_1=NULL
head(full_heat_tablemelt1)
heat2 = dcast(full_heat_tablemelt, gene_1 ~ gP1, fun.aggregate = mean, value.var = 'pvalue')
row.names(heat2) = heat2$gene_1
heat2$gene_1=NULL
heat2anno = ifelse(heat2<0.05, '*', '')

anno_table =as.data.frame(row.names(heat1))
row.names(anno_table) = anno_table$`row.names(heat1)`
rownames(anno_table) = anno_table$`row.names(heat1)`
anno_table$pathway = full_paths1$pathways[match(anno_table$`row.names(heat1)`, full_paths1$gene_symbol)]

anno_table1 = anno_table 
anno_table1$`row.names(heat1)`=NULL
pdf(file = 'heatmap of pathway and gene correlation in each cell type - full.pdf')
pheatmap::pheatmap(heat1, display_numbers = heat2anno, annotation_row = anno_table1, fontsize_row = 5, color = colorRampPalette(c( "#9a133d", "#ffffff", '#1a318b'))(100))
dev.off()





