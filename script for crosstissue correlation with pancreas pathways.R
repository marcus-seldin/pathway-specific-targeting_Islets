setwd('G:/My Drive/lab files/pancreas targeting/GTEx islet enrichments/datasets')
sex_table = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')


#install.packages('BiocManager')
library(WGCNA)
#BiocManager::install('WGCNA')
library(ggplot2)
install.packages('Rmisc')

library(hexbin)
library(dplyr)
library(ggrepel)
library(Rmisc)
library(reshape2)
allowWGCNAThreads()
#Try top percentile cor across tissues
load('G:/My Drive/lab files/pancreas targeting/GTEx islet enrichments/datasets/GTEx NA included env.RData')

working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
working_dataset[1:5,1:5]
sex_table = read.delim('G:/My Drive/lab files/pancreas targeting/GTEx islet enrichments/datasets/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
table(sex_table$sexMF)

new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
table(new_trts$sexMF)


#read in Secreted Proteins
Secreted_proteins <- read.delim("G:/My Drive/lab files/pancreas targeting/GTEx islet enrichments/datasets/uniprot-secreted-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab", header = T, check.names = F)
Secreted_proteins[1:10,1:5]

deg_table = read.delim('G:/My Drive/lab files/pancreas targeting/GTEx islet enrichments/datasets/uniprot-human-genes and goterms mapping.tab')

path_set = deg_table[grepl('insulin secretion', deg_table$Gene.ontology..biological.process.),]

tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]


tissue2 <- working_dataset[,grepl('Pancreas', colnames(working_dataset)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue2 = tissue2[,colnames(tissue2) %in% path_set$Gene.names...primary..]
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')

tt1 = tissue.tissue.p$p
tt1[is.na(tt1)] = 0.5
tt1[tt1==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))

colnames(cc3) = 'Ssec_score'
cc3$gene_tissue = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)
scores=cc3
mean(scores$Ssec_score + (sd(scores$Ssec_score)*3))
scores$gene_symbol =  gsub("\\_.*","",scores$gene_tissue)
scores$secreted = ifelse(scores$gene_symbol %in% Secreted_proteins$`Gene names  (primary )`, 'Secreted', 'Non-secreted')
table(scores$secreted)
scores$Sig_P1e3 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*3)), 'Significant', 'NS')
tissue_col_map = as.data.frame(table(gsub(".*_","",scores$gene_tissue)))
tissue_col_map$col = c('darkorange2', 'blue1', 'firebrick1', 'darkorchid1')

scores$tissue_col = tissue_col_map$col[match(gsub(".*_","",scores$gene_tissue), tissue_col_map$Var1)]
scores$sec_color = ifelse(scores$secreted=='Secreted', 'darkorange1', 'darkorchid')

scores1 = scores[scores$Sig_P1e3=='Significant',]
table(gsub(".*_","",scores1$gene_tissue))

pie(table(gsub(".*_","",scores1$gene_tissue)), main="Tissue distribution of significant cors", cex=0.8)

scores1 = scores1[order(scores1$Ssec_score, decreasing = T),]
scores2 = scores1[1:100,]
pie(table(gsub(".*_","",scores2$gene_tissue)), main="Tissue distribution of top 100")

scores2 = scores1[1:50,]
library(forcats)
ggplot(scores2, aes(x=fct_reorder2(gene_tissue, Ssec_score, Ssec_score, .desc = T), y=Ssec_score)) + geom_col(fill=scores2$sec_color) +  ggtitle('Top genes enriched for Insulin secretion pathways') + theme_minimal()  + xlab('gene') + ylab('cross-tissue score (Ssec) in GTEx')+ theme(axis.text.x = element_text(angle=90, size=8,vjust =0.5, hjust = 0.5), plot.title = element_text(hjust=0.5)) 

ss1 = scores[scores$secreted=='Secreted',]
scores2 = ss1[1:40,]

ggplot(scores2, aes(x=fct_reorder2(gene_tissue, Ssec_score, Ssec_score, .desc = T), y=Ssec_score)) + geom_col(fill=scores2$sec_color) +  ggtitle('Top SECRETED genes enriched for Insulin secretion pathways') + theme_minimal()  + xlab('gene') + ylab('cross-tissue score (Ssec) in GTEx')+ theme(axis.text.x = element_text(angle=90, size=8,vjust =0.5, hjust = 0.5), plot.title = element_text(hjust=0.5)) 





#look at relative enrichemnt within vs across tissues
top_gene = ss1$gene_symbol[4]
top_gene_target = paste0(top_gene, '_', 'Pancreas')
top_gene_ttissues = 'Pancreas'
top_gene_origin = ss1$gene_tissue[4]
top_gene_otissues = gsub(".*_","",ss1$gene_tissue[4])


#this should not require editing to run
fc_table = scores
fc_table$tissues = gsub(".*_","",fc_table$gene_tissue)
fc_table1 = fc_table[fc_table$tissues==top_gene_otissues,]
fc_table1$normSsec = fc_table1$Ssec_score/mean(fc_table1$Ssec_score)

binned_com_table = fc_table1 %>% dplyr::select(gene_tissue, normSsec)
binned_com_table$tissue_cat = paste0(top_gene_otissues, '-', top_gene_ttissues)

#look at within target tissue comparisons
tissue1 <- working_dataset[,grepl(top_gene_ttissues, colnames(working_dataset)),]
tissue2 <- working_dataset[,grepl(top_gene_ttissues, colnames(working_dataset)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue2 = tissue2[,colnames(tissue2) %in% path_set$Gene.names...primary..]
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')

tt1 = tissue.tissue.p$p
tt1[is.na(tt1)] = 0.5
tt1[tt1==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))

colnames(cc3) = 'Ssec_score'
cc3$gene_tissue = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
cc3$normSsec = cc3$Ssec_score/mean(cc3$Ssec_score)
new_table = cc3 %>% dplyr::select(gene_tissue, normSsec)
new_table$tissue_cat = paste0(top_gene_ttissues, '-', top_gene_ttissues)
binned_com_table = as.data.frame(rbind(binned_com_table, new_table))


tissue1 <- working_dataset[,grepl(top_gene_otissues, colnames(working_dataset)),]
tissue2 <- working_dataset[,grepl(top_gene_otissues, colnames(working_dataset)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue2 = tissue2[,colnames(tissue2) %in% path_set$Gene.names...primary..]
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')

tt1 = tissue.tissue.p$p
tt1[is.na(tt1)] = 0.5
tt1[tt1==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))

colnames(cc3) = 'Ssec_score'
cc3$gene_tissue = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
cc3$normSsec = cc3$Ssec_score/mean(cc3$Ssec_score)
new_table = cc3 %>% dplyr::select(gene_tissue, normSsec)
new_table$tissue_cat = paste0(top_gene_otissues, '-', top_gene_otissues)
binned_com_table = as.data.frame(rbind(binned_com_table, new_table))

binned_com_table$gene_symbol =  gsub("\\_.*","",binned_com_table$gene_tissue)


#now compare
sec_comp = binned_com_table[binned_com_table$gene_symbol==top_gene,]

ggplot(sec_comp, aes(x=tissue_cat, y=normSsec, fill=tissue_cat)) + geom_col(show.legend = FALSE) + theme_classic() + scale_fill_manual(values = c('brown1', 'darkorchid2', 'dodgerblue3')) + ylab('Normalized Ssec score') + xlab('') + ggtitle(paste0('Within vs across scores ', top_gene, ' ', top_gene_otissues, '-', top_gene_ttissues))
