# Code is used from http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# and from https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

# DESeq2 of nondiscovery
library(DESeq2)
library(ggplot2)
# load data ####
setwd("C:/Users/WSY/Desktop/My_data/From_HPC/raw_count_matrix")
count_gene_nondis = read.csv("gene_count_matrix_nondiscovery_10q.csv",row.names = 1, header = T)
colTable = read.csv("experiment_design_table.csv", stringsAsFactors = T, row.names = 1, header = T)

# create DESeq2 object ####
dds = DESeqDataSetFromMatrix(countData = count_gene_nondis, colData = colTable, design = ~ Condition)
dds
# 55414 genes of nondiscovery in total

# pre-filter the gene list ####
valid_dds = (rowSums(counts(dds))>0)
dds = dds[valid_dds,]
dim(dds) 
#27778

# Exploratory analysis and visualization ####
# variance stabilizing transformation (VST) and rlog
rld = rlog(dds, blind = F)
vsd = vst(dds, blind = F)
dds_test = estimateSizeFactors(dds)
sizeFactors(dds_test)
normalized_counts = counts(dds_test, normalized=TRUE)
write.csv(normalized_counts, "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/gene_counts.csv")

# boxplot-to compare between data with and without normalization
lgc_raw_dds = log2(counts(dds, normalized=F) + 1)
lgc_norm_dds = log2(normalized_counts + 1)
par(mfrow = c(2,1))
boxplot(lgc_raw_dds, main = "log2(raw count of nondiscovery)")
boxplot(lgc_norm_dds, main = "log2(normalized count of nondiscovery)")

# plot densities
library("limma")
par(mfrow = c(1,1))
plotDensities(lgc_raw_dds, group = dds$Condition, legend = "topright")
plotDensities(lgc_norm_dds, group = dds$Condition, legend = "topright")
# slightly different 

# PCA
plotPCA(rld, intgroup="Condition", returnData=F) +
  geom_text(aes(label=rld@colData@listData[["Sample"]]), vjust=-0.6, hjust=1.1, fontface=0.5, size=3)

plotPCA(vsd, intgroup="Condition") +
  geom_text(aes(label=rld@colData@listData[["Sample"]]), vjust=-0.5, hjust=0.8, fontface=0.5, size=3)

# distance plot
distance_dds = dist(t(assay(vsd)))
library("RColorBrewer")
library("DESeq2")
library("pheatmap")
DistMatrix = as.matrix(distance_dds)
rownames(DistMatrix) = paste(vsd$Condition, vsd$Sample.1, sep="-")
colnames(DistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(DistMatrix,
         clustering_distance_rows=distance_dds,
         clustering_distance_cols=distance_dds,
         col=colors)

# do DE analysis ####
dds = DESeq(dds)

# dispersion plot ####
plotDispEsts(dds)

# build the result tables ####
res_KO = results(dds, contrast = c("Condition", "KOHS", "KONS"))
summary(res_KO)
res_KO = res_KO[order(res_KO$padj),]
head(res_KO)
res_sig_KO = subset(res_KO,res_KO$padj < 0.05 & abs(res_KO$log2FoldChange) > 1)
FDR_KO = subset(res_KO,res_KO$padj < 0.05)
write.csv(as.data.frame(res_sig_KO), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_KOHSvsKONS_nondiscovery_10q.csv")
write.csv(as.data.frame(res_KO), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KOHSvsKONS_nondiscovery_10q.csv")
write.csv(as.data.frame(FDR_KO), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/FDR_KOHSvsKONS_nondiscovery_10q.csv")

res_Wt = results(dds, contrast = c("Condition", "WtHS", "WtNS"))
summary(res_Wt)
res_Wt = res_Wt[order(res_Wt$padj),]
head(res_Wt)
res_sig_Wt = subset(res_Wt,res_Wt$padj < 0.05 & abs(res_Wt$log2FoldChange) > 1)
FDR_Wt = subset(res_Wt,res_KO$padj < 0.05)
write.csv(as.data.frame(res_sig_Wt), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_WtHSvsWtNS_nondiscovery_10q.csv")
write.csv(as.data.frame(res_Wt), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/WtHSvsWtNS_nondiscovery_10q.csv")
write.csv(as.data.frame(FDR_Wt), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/FDR_WtHSvsWtNS_nondiscovery_10q.csv")

res_HS = results(dds, contrast = c("Condition", "KOHS", "WtHS"))
summary(res_HS)
res_HS = res_HS[order(res_HS$padj),]
head(res_HS)
FDR_HS = subset(res_HS,res_HS$padj < 0.05)
res_sig_HS = subset(res_HS,res_HS$padj < 0.05 & abs(res_HS$log2FoldChange) > 1)
write.csv(as.data.frame(res_sig_HS), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_KOHSvsWtHS_nondiscovery_10q.csv")
write.csv(as.data.frame(res_HS), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KOHSvsWtHS_nondiscovery_10q.csv")
write.csv(as.data.frame(FDR_HS), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/FDR_KOHSvsWtHS_nondiscovery_10q.csv")

res_NS = results(dds, contrast = c("Condition", "KONS", "WtNS"))
summary(res_NS)
res_NS = res_NS[order(res_NS$padj),]
head(res_NS)
FDR_NS = subset(res_NS,res_NS$padj < 0.05)
res_sig_NS = subset(res_NS,res_NS$padj < 0.05 & abs(res_NS$log2FoldChange) > 1)
write.csv(as.data.frame(res_sig_NS), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_KONSvsWtNS_nondiscovery_10q.csv")
write.csv(as.data.frame(res_NS), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KONSvsWtNS_nondiscovery_10q.csv")
write.csv(as.data.frame(FDR_NS), file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/FDR_KONSvsWtNS_nondiscovery_10q.csv")


# annotation result table ####
library("AnnotationDbi")
library(biomaRt)
library('org.Mm.eg.db')

ensembl=useMart('ensembl')
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
attributes = listAttributes(ensembl)
filters = listFilters(ensembl)

file_list = list("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KOHSvsKONS_nondiscovery_10q.csv",
                 "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/WtHSvsWtNS_nondiscovery_10q.csv",
                 "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KOHSvsWtHS_nondiscovery_10q.csv",
                 "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KONSvsWtNS_nondiscovery_10q.csv")

for (file in file_list) {
  print(file)
  data = read.csv(file = file)
  print(head(data))
  gene_id = data$X
  info =  getBM(attributes = c("ensembl_gene_id","external_gene_name", "chromosome_name", 
                               "start_position", "end_position", "gene_biotype", "entrezgene_id", "description"),
                filters = "ensembl_gene_id",
                values = gene_id,
                mart = ensembl)
  merge = merge(data, info, by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE, all.y = T)
  merge = merge[order(merge$padj),]
  write.csv(merge, file = file)
} 


# plotting result
# MA-plot
pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KO_MA.pdf",  
    width = 8, height = 8)
DESeq2::plotMA(res_KO, alpha=0.05, xlab = "mean of normalized counts of KOHSvsKONS in nondiscovery")
abline(h=c(-1,1), col="black")
dev.off()

pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/Wt_MA.pdf",  
    width = 8, height = 8)
DESeq2::plotMA(res_Wt, alpha=0.05, xlab = "mean of normalized counts of WtHSvsWtNS in nondiscovery")
abline(h=c(-1,1), col="black")
dev.off()

pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/HS_MA.pdf",  
    width = 8, height = 8)
DESeq2::plotMA(res_HS, alpha=0.05, xlab = "mean of normalized counts of KOHSvsWtHS in nondiscovery")
abline(h=c(-1,1), col="black")
dev.off()

pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/NS_MA.pdf",  
    width = 8, height = 8)
DESeq2::plotMA(res_NS, alpha=0.05, xlab = "mean of normalized counts of KONSvsWtNS in nondiscovery")
abline(h=c(-1,1), col="black")
dev.off()

# count plots
# UMOD = ENSMUSG00000030963
plotCounts(dds, gene = 'ENSMUSG00000030963', intgroup=c("Condition"))
library("ggbeeswarm")

geneCounts <- plotCounts(dds, gene = 'ENSMUSG00000030963', intgroup = c("Condition","Sample"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Condition, y = count, color = Sample)) +  geom_beeswarm(cex = 3) + labs(title = "UMOD level") + 
  geom_text(aes(label=dds@colData@listData[["Sample"]]), vjust=-0.5, hjust=0.8, fontface=0.5, size=3)


# volcano plot ####
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

HS_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/res_with_annotation/KOHSvsWtHS_nondiscovery_10q.csv")
pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/HS_volcano.pdf",  
    width = 10, height = 8)
EnhancedVolcano(HS_file, lab = HS_file$external_gene_name, 
                x = 'log2FoldChange', y = 'padj', title = "Non-Discovery-10_KOHS_vs_WtHS",
                ylab = bquote(~-Log[10] ~ italic(P-adj)),
                axisLabSize = 15,
                titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 6,legendLabSize = 10,
                legendIconSize = 5, pCutoff = 5e-2, pCutoffCol = 'padj',
                FCcutoff = 1, pointSize = 2, labSize = 5, 
                drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, 
                legendLabels=c('Not sig.',expression(Log[2] ~ FC),'p-adj',
                               expression(p-adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'top')
dev.off()

KO_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/res_with_annotation/KOHSvsKONS_nondiscovery_10q.csv")
pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/KO_volcano.pdf",  
    width = 10, height = 8)
EnhancedVolcano(KO_file, lab = KO_file$external_gene_name, 
                x = 'log2FoldChange', y = 'padj', title = "Non-Discovery-10_KOHS_vs_KONS",
                ylab = bquote(~-Log[10] ~ italic(P-adj)),
                axisLabSize = 15,
                titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 6,legendLabSize = 10,
                legendIconSize = 5, pCutoff = 5e-2, pCutoffCol = 'padj',
                FCcutoff = 1, pointSize = 2, labSize = 5, 
                drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, 
                legendLabels=c('Not sig.',expression(Log[2] ~ FC),'p-adj',
                               expression(p-adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'top')
dev.off()


Wt_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/res_with_annotation/WtHSvsWtNS_nondiscovery_10q.csv")
pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/Wt_volcano.pdf",  
    width = 10, height = 8)
EnhancedVolcano(Wt_file, lab = Wt_file$external_gene_name, 
                x = 'log2FoldChange', y = 'padj', title = "Non-Discovery-10_WtHS_vs_WtHS",
                ylab = bquote(~-Log[10] ~ italic(P-adj)),
                axisLabSize = 15,
                titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 6,legendLabSize = 10,
                legendIconSize = 5, pCutoff = 5e-2, pCutoffCol = 'padj',
                FCcutoff = 1, pointSize = 2, labSize = 5, 
                drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, 
                legendLabels=c('Not sig.',expression(Log[2] ~ FC),'p-adj',
                               expression(p-adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'top')
dev.off()

NS_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/res_with_annotation/KONSvsWtNS_nondiscovery_10q.csv")
pdf(file = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/NS_volcano.pdf",  
    width = 10, height = 8)
EnhancedVolcano(NS_file, lab = NS_file$external_gene_name, 
                x = 'log2FoldChange', y = 'padj', title = "Non-Discovery-10_KONS_vs_WtNS",
                ylab = bquote(~-Log[10] ~ italic(P-adj)),
                axisLabSize = 15,
                titleLabSize = 10, subtitleLabSize = 8,
                captionLabSize = 6,legendLabSize = 10,
                legendIconSize = 5, pCutoff = 5e-2, pCutoffCol = 'padj',
                FCcutoff = 1, pointSize = 2, labSize = 5, 
                drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 50, 
                legendLabels=c('Not sig.',expression(Log[2] ~ FC),'p-adj',
                               expression(p-adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'top')
dev.off()

# volcano plot with transport genes ####
# generate a table only contains transporter genes
colnames(HS_file)[1] = "HS"
colnames(NS_file)[1] = "NS"
colnames(KO_file)[1] = "KO"
colnames(Wt_file)[1] = "Wt"
file_list = list(HS_file, NS_file, KO_file, Wt_file)

library(dplyr)
for (ele in file_list){
  transporter =dplyr::filter(ele, grepl('sodium|channel|transporter|solute|transmembrane', description))
  transporter$type = 'transporter'
  non_transpoter = ele[!ele$X %in% transporter$X,]
  non_transpoter$type = "non-transporter"
  merge_file = rbind(transporter, non_transpoter)
  merge_file = merge_file %>%
    mutate(sig_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                log2FoldChange <= -1 & padj <= 0.05 ~ "down", TRUE ~ "ns"))
  sig_table =dplyr::filter(merge_file, grepl('down|up', sig_type))
  
  if (colnames(ele)[1] == "HS"){
    new_HS = merge_file
    transporter_HS = transporter
    sig_HS = sig_table
  } else if (colnames(ele)[1] == "NS") {
    new_NS = merge_file
    transporter_NS = transporter
    sig_NS = sig_table
  } else if (colnames(ele)[1] == "KO") {
    new_KO = merge_file
    transporter_KO = transporter
    sig_KO = sig_table
  } else {
    new_Wt = merge_file
    transporter_Wt = transporter
    sig_Wt = sig_table
  }
} 

cols = c("up" = "#0077b6", "down" = "#ae2012", "not sig" = "#999999") 
#alphas = c("up" = 0.6, "down" = 0.6, "ns" = 0.5)
sizes = c("up" = 2.5, "down" = 2.5, "not sig" = 2)

# plot - 1 
vol_plot_1 =  ggplot(data = new_HS, aes(x = log2FoldChange,
             y = -log10(padj), fill = sig_type, size = sig_type)) + 
            geom_point(colour = "black", shape = 21, alpha=0.6) + 
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
            geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
            scale_fill_manual(values = cols)+
            scale_size_manual(values = sizes)+
            geom_point(data = transporter_HS,       
                        size = 3, shape = 21,fill = "#fcbf49", colour = "black")
vol_plot_1

# plot - 2
vol_plot_2 =  ggplot(data = new_Wt, aes(x = log2FoldChange, # change name
                                      y = -log10(padj))) + 
  geom_point(aes(fill = sig_type), 
             alpha = 0.7, shape = 21, size = 3) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_label_repel(data = sig_Wt, # change name
                   aes(label = external_gene_name),
                   force = 1,
                   nudge_y = 0.2,max.overlaps = 25)+
  labs(title = "Gene expression changes in WtNS vs KONS samples",# change name
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

vol_plot_2

# plot - 3 only label transporter genes
transport_vol_plot =  ggplot(data = new_Wt, aes(x = log2FoldChange, # change name
                                      y = -log10(padj))) + 
  geom_point(aes(fill = sig_type), 
             alpha = 0.3, shape = 21, size = 2, fill = "#999999") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_point(data = transporter_Wt,       # change name
             size = 3, shape = 21, alpha = 1.5 ,fill = "#fcbf49", colour = "black")+
  geom_label_repel(data = transporter_Wt, # change name
                   aes(label = external_gene_name),
                   force = 1,
                   nudge_y = 0.5, max.overlaps = 50)+
  labs(title = "Gene expression changes in WtHS vs KOHS samples",# change name
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

transport_vol_plot

# gene clustering
library("genefilter")
library(DESeq2)
library(pheatmap)
topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  = assay(vsd)[topVarGenes, ]
mat  = mat - rowMeans(mat)
pheatmap(mat, annotation_col = colTable)

# make gene expression level plot for gene with adjp=1 ####
library("ggbeeswarm")
adjp_1_HS = filter(HS_file, padj==1)
adj1_HS_name = c(adjp_1_HS$X)
adj1_HS_symbol = c(adjp_1_HS$external_gene_name)

for (gene in adj1_HS_name[1:50]){
  id = which(adj1_HS_name == gene)
  symbol = adj1_HS_symbol[id]
  pdf(paste("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/gene_level/", 
            symbol, ".pdf", sep = ""), width = 6, height = 5)
  geneCounts = plotCounts(dds, gene = gene, intgroup = c("Condition","Sample.1"),returnData = TRUE)
  plot = ggplot(geneCounts, aes(x = Condition, y = count, color = Sample.1)) +
    scale_y_log10() +  geom_beeswarm(cex = 3) + labs(title = symbol) + 
    geom_text(label=colTable$Sample.1, 
              nudge_x = 0.15, nudge_y = 0, 
              check_overlap = T)
  pdf("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/gene_level/", 
      paste(symbol, ".pdf", sep = ""), width = 6, height = 5)
  print(plot)
  dev.off()
}

for (ele in range(1,50)) {
  dev.off()
}
