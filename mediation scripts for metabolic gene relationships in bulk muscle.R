setwd('G:/My Drive/lab files/Lauren Albrecht/muslce pfk story/mediation examples')
load('G:/My Drive/lab files/sex-difference myokine study/github/raw files/GTEx NA included env.RData')
human_pathways = read.delim('G:/My Drive/Datasets/Human/genome files/uniprot-human-genes and goterms mapping.tab')
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
filtered_working = GTEx_subfiltered
filtered_working$gene_tissue=NULL

lauren_paths = read.csv('G:/My Drive/lab files/Lauren Albrecht/muslce pfk story/lauren curated pathways.csv')
table(lauren_paths$Pathway)
table(lauren_paths$Gene[lauren_paths$Pathway=='Lysosome'])


filtered_working = GTEx_full
table(colnames(filtered_working[grepl('ITIH5', colnames(filtered_working))]))

gene_tissue_origin = 'PFKM_Muscle - Skeletal'
covariate_gene_tissue = 'LAMP1_Muscle - Skeletal'
pathway_term = 'PPP'
pathway_tissue = 'Muscle - Skeletal'


run_full_adj = function(gene_tissue_origin, covariate_gene_tissue, pathway_term, pathway_tissue){
  
  new_paths = lauren_paths[lauren_paths$Pathway==pathway_term,]
  new_paths$tissuepaste = paste0(pathway_tissue)
  new_paths$gene_tissue = paste0(new_paths$Gene, '_', new_paths$tissuepaste)
  #write.csv(new_paths, file = paste0('Pathway genes listed for ', pathway_tissue, ' ', pathway_term, '.csv'), row.names = F)
  
  gg1 = as.data.frame(filtered_working[,colnames(filtered_working) %in% gene_tissue_origin])
  row.names(gg1) = row.names(filtered_working)
  gg2 = as.data.frame(filtered_working[,colnames(filtered_working) %in% new_paths$gene_tissue])
  row.names(gg2) = row.names(filtered_working)
  gg3 = as.data.frame(filtered_working[,colnames(filtered_working) %in% covariate_gene_tissue])
  row.names(gg3) = row.names(filtered_working)
  
  combined_key = as.data.frame(cbind(gg1, gg3))
  combined_key = na.omit(combined_key)
  gg2 = gg2[row.names(gg2) %in% row.names(combined_key),]
  
  cor1 = bicorAndPvalue(combined_key$`filtered_working[, colnames(filtered_working) %in% gene_tissue_origin]`, gg2)
  resids = summary(lm(combined_key$`filtered_working[, colnames(filtered_working) %in% gene_tissue_origin]` ~ combined_key$`filtered_working[, colnames(filtered_working) %in% covariate_gene_tissue]`))$residuals
  
  
  cor2 = bicorAndPvalue(resids, gg2)
  heat_set = as.data.frame(rbind(cor1$bicor, cor2$bicor))
  heat_setp = as.data.frame(rbind(cor1$p, cor2$p))
  row.names(heat_set) = c(paste0(gene_tissue_origin, '-cor'), paste0(covariate_gene_tissue, ' ADJ-cor'))
  row.names(heat_setp) = c(paste0(gene_tissue_origin, '-cor'), paste0(covariate_gene_tissue, '-adj-cor'))
  row.names(heat_set) = gsub(paste0('_', pathway_tissue), '', row.names(heat_set))
  colnames(heat_set) = gsub(paste0('_', pathway_tissue), '', colnames(heat_set))
  pval_settings = ifelse(heat_setp<0.01, '*', '')
  pval_settings = ifelse(heat_setp<1e-3, '**', pval_settings)
  pval_settings = ifelse(heat_setp<1e-6, '***', pval_settings)
  fontsize_number = 20 ## increase label size
  pdf(file = paste0('heatmap of ', gene_tissue_origin, ' correlated with ', pathway_tissue, ' ', pathway_term, ' ', covariate_gene_tissue, '-adjusted.pdf'))
  pp1 = pheatmap(heat_set, cluster_rows = F, cluster_cols = F, display_numbers = pval_settings, fontsize_number = 30, color = colorRampPalette(c("darkorchid", "darkorange", "white"))(100))
  print(pp1)
  dev.off()
  
  pmelted =melt(as.matrix(heat_setp))
  head(pmelted)
  colnames(pmelted) = c('normal_adj', 'gene_tissue', 'pvalue')
  pmelted$log10p = -log10(pmelted$pvalue)
  
  pdf(file = paste0('boxplot comparison of ', gene_tissue_origin, ' correlated with ', pathway_tissue, ' ', pathway_term, ' ', covariate_gene_tissue, '-adjusted.pdf'))
  pp2 = ggboxplot(pmelted, x = "normal_adj", y = "log10p", 
                  color = "normal_adj", palette = c("darkorange", "darkorchid")
  ) + stat_compare_means() + geom_point(aes(x=normal_adj, y=log10p, fill = normal_adj, col=normal_adj), position = position_jitterdodge(0.1, dodge.width = 0.8)) 
  print(pp2)
  dev.off()
  head(pmelted)
  
  pmelted =melt(as.matrix(heat_set))
  head(pmelted)
  colnames(pmelted) = c('normal_adj', 'gene_tissue', 'cor_coeff')
  
  
  pdf(file = paste0('Directional boxplot comparison of ', gene_tissue_origin, ' correlated with ', pathway_tissue, ' ', pathway_term, ' ', covariate_gene_tissue, '-adjusted.pdf'))
  pp2 = ggboxplot(pmelted, x = "normal_adj", y = "cor_coeff", 
                  color = "normal_adj", palette = c("darkorange", "darkorchid")
  ) + stat_compare_means() + geom_point(aes(x=normal_adj, y=cor_coeff, fill = normal_adj, col=normal_adj), position = position_jitterdodge(0.1, dodge.width = 0.8)) 
  print(pp2)
  dev.off()
}



run_full_adj('PFKM_Muscle - Skeletal', 'LAMP1_Muscle - Skeletal', 'PPP', 'Muscle - Skeletal')
run_full_adj('PFKM_Muscle - Skeletal', 'LAMP2_Muscle - Skeletal', 'PPP', 'Muscle - Skeletal')
