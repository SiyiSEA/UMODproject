# Part of code is used from https://github.com/gencorefacility/r-notebooks/blob/master/gsea.Rmd

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)

# set the desired organism here
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

# prepare input data ####
# reading in data
# load genes files in each contrast
setwd("C:/Users/WSY/Desktop/My_data/GSEA")
GO_HS_file = read.csv("./pre_data/KOHSvsWtHS_nondiscovery_10q.csv", header = T)
colnames(GO_HS_file)[1] = "HS_gene_id"
FDR_HS_file = subset(GO_HS_file, GO_HS_file$padj < 0.05)
SIG_HS_file = subset(GO_HS_file, GO_HS_file$padj < 0.05 & abs(GO_HS_file$log2FoldChange) > 1)

GO_NS_file = read.csv("./pre_data/KONSvsWtNS_nondiscovery_10q.csv", header = T)
colnames(GO_NS_file)[1] = "NS_gene_id"
FDR_NS_file = subset(GO_NS_file,GO_NS_file$padj < 0.05)
SIG_NS_file = subset(GO_NS_file,GO_NS_file$padj < 0.05 & abs(GO_NS_file$log2FoldChange) > 1)

GO_Wt_file = read.csv("./pre_data/WtHSvsWtNS_nondiscovery_10q.csv", header = T)
colnames(GO_Wt_file)[1] = "Wt_gene_id"
FDR_Wt_file = subset(GO_Wt_file,GO_Wt_file$padj < 0.05)
SIG_Wt_file = subset(GO_Wt_file,GO_Wt_file$padj < 0.05 & abs(GO_Wt_file$log2FoldChange) > 1)

GO_KO_file = read.csv("./pre_data/KOHSvsKONS_nondiscovery_10q.csv", header = T)
colnames(GO_KO_file)[1] = "KO_gene_id"
FDR_KO_file = subset(GO_KO_file,GO_KO_file$padj < 0.05)
SIG_KO_file = subset(GO_KO_file,GO_KO_file$padj < 0.05 & abs(GO_KO_file$log2FoldChange) > 1)


# GO ####
# for GO 
all_genes_list = list(GO_HS_file, GO_NS_file, GO_Wt_file, GO_KO_file)
setwd("C:/Users/WSY/Desktop/My_data/GSEA/GO_Pathways_all_genes/")

FDR_genes_list = list(FDR_HS_file, FDR_NS_file, FDR_Wt_file, FDR_KO_file)
setwd("C:/Users/WSY/Desktop/My_data/GSEA/GO_Pathway_FDR_genes/")

SIG_genes_list = list(SIG_HS_file, SIG_NS_file, SIG_Wt_file, SIG_KO_file)
setwd("C:/Users/WSY/Desktop/My_data/GSEA/GO_Pathway_sig_genes/")

# check the first colname
for (df in FDR_genes_list) {
  print(colnames(df)[1])
}

# change file list
for (df in FDR_genes_list) {
  if (colnames(df)[1] == "HS_gene_id") {
    path = "./HS/"
    name = "HS"
    title_name = "KOHS vs WtHS Enriched Pathways GO"
    GO_filename = "./HS/HS_GO.tsv"
    GO_object = "./HS/HS_GO.rds"
    plot_size = c(6,3,6,3,6,3,8,5,6,6,6,6,6,5,6,3,6,6)
  } else if (colnames(df)[1] == "NS_gene_id") {
    path = "./NS/"
    name = "NS"
    title_name = "KONS vs WtNS Enriched Pathways GO"
    GO_filename = "./NS/NS_GO.tsv"
    GO_object = "./NS/NS_GO.rds"
    plot_size = c(7,7,6,4,8,8,13,15,6,6,14,10,20,16,10,12,6,6)
  } else if (colnames(df)[1] == "Wt_gene_id") {
    path = "./Wt/"
    name = "Wt"
    title_name = "WtHS vs WtNS Enriched Pathways GO"
    GO_filename = "./Wt/Wt_GO.tsv"
    GO_object = "./Wt/Wt_GO.rds"
    plot_size = c(8,6,6,4,8,8,16,5,6,6,5,5,5,4,10,6,6,6)
  } else {
    path = "./KO/"
    name = "KO"
    title_name = "KOHS vs KONS Enriched Pathways GO"
    GO_filename = "./KO/KO_GO.tsv"
    GO_object = "./KO/KO_GO.rds"
    plot_size = c(8,8,6,4,8,8,16,5,8,8,10,10,20,18,10,10,6,6)
  }
  
  print(path)
  print(name)
  print(title_name)
  print(GO_filename)

  # prepare data
  original_gene_list = df$log2FoldChange
  class(original_gene_list)
  names(original_gene_list) = df[, c(1)]
  class(names(original_gene_list))
  # View(original_gene_list)
  gene_list = na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = T)
  
  # check which options are available 
  keytypes(org.Mm.eg.db)
  gse = gseGO(geneList = gene_list, ont = "ALL", keyType = "ENSEMBL", nPerm = 10000,
              minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, 
              verbose = T, OrgDb = org.Mm.eg.db, pAdjustMethod = "none")
  
  # save gse object
  saveRDS(gse, file = GO_object)
  
  # add -log(pvalue) as a column
  gse@result$"-log(p-value)" = -log(gse@result$pvalue)
  gse@result = gse@result[order(-gse@result$`-log(p-value)`),]
  
  gse_PC = setReadable(gse, OrgDb = organism) # some gene id do not have a symbol??
  #write.table(gse, file=GO_filename, sep="\t", quote=F, row.names = F)
  
  # dotplot
  require(DOSE)
  dotplot_1 = dotplot(gse, showCategory=10, title = title_name, 
                      split=".sign",label_format=30, font.size = 10) + facet_grid(.~.sign)
  pdf(paste(path, "dotplot_1_",name, ".pdf", sep = ""), width = plot_size[1], height = plot_size[2])
  print(dotplot_1)
  dev.off()
  
  dotplot_2 = dotplot(gse, showCategory=10, title = title_name, label_format=30, font.size = 10) 
  pdf(paste(path, "dotplot_2_",name, ".pdf", sep = ""), width = plot_size[3], height = plot_size[4])
  print(dotplot_2)
  dev.off()
  
  dotplot_3 = dotplot(gse, showCategory=10, title = title_name, 
                      split="ONTOLOGY",label_format=30, font.size = 10) +
    facet_grid(ONTOLOGY~., scale='free')+
    scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
  pdf(paste(path, "dotplot_3_",name, ".pdf", sep = ""), width = plot_size[5], height = plot_size[6])
  print(dotplot_3)
  dev.off()
  
  # barplot
  ggdata = data.frame(value = gse@result$`-log(p-value)`[1:10],
                      pathway = gse@result$Description[1:10],
                      ontology = gse@result$ONTOLOGY[1:10])
  
  barplot = ggplot(ggdata, aes(x = reorder(pathway, +value), y = value,  fill = ontology)) +    
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    labs(x = "Pathway", y = "-log(p-value)", title = title_name) +
    theme(axis.text.x = element_text(angle = 90, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10)) +
    geom_text(
      aes(label = paste(round(value, 3))),
      color = "black",size = 4,hjust=-0.1,
      position = position_dodge(0.9)) + coord_flip() +
    scale_fill_manual(values=c("#5fad56", "#f2c14e", "#f78154"))
  
  pdf(paste(path, "barplot_", name, ".pdf", sep = ""), width = plot_size[7], height = plot_size[8])
  print(barplot)
  dev.off()

  # enrichment map
  x2 = pairwise_termsim(gse)
  enrichment_plot = emapplot(x2, showCategory = 10)
  pdf(paste(path, "emapplot_", name, ".pdf", sep = ""), width = plot_size[9], height = plot_size[10])
  print(enrichment_plot)
  dev.off()
  
  # category netplot
  # categorySize can be either 'pvalue' or 'geneNum'
  netplot = cnetplot(gse_PC, categorySize="pvalue", foldChange=gene_list, 
           showCategory = 3, cex_label_gene = 0.5)
  pdf(paste(path, "netplot_", name, ".pdf", sep = ""), width = plot_size[11], height = plot_size[12])
  print(netplot)
  dev.off()
  
  # circular centplot
  circular_cen = cnetplot(gse_PC, categorySize="pvalue", foldChange=gene_list, 
           showCategory = 3, cex_label_gene = 0.5, circular = T, colorEdge = TRUE)
  pdf(paste(path, "circular_cenplot_",name, ".pdf", sep = ""), width = plot_size[13], height = plot_size[14])
  print(circular_cen)
  dev.off()
  
  
  # heatmap
  heatplot(gse_PC, showCategory = 10, foldChange=gene_list, label_format = 30)
  
  # Ridgeplot
  library(ggridges)
  ridge = ridgeplot(gse, showCategory = 10,
            fill = "p.adjust",label_format = 20) + labs(x = "enrichment distribution")
  pdf(paste(path, "ridgeplot_", name, ".pdf", sep = ""), width = plot_size[15], height = plot_size[16])
  print(ridge)
  dev.off()
  
  #GSEA Plot
  gsea = gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
  pdf(paste(path, "gseaplot_", name, ".pdf", sep = ""), width = plot_size[17], height = plot_size[18])
  print(gsea)
  dev.off()

}


# KEGG ####

# for KEGG (choose one of lists and corresponding path)
all_genes_list = list(GO_HS_file, GO_NS_file, GO_Wt_file, GO_KO_file)
setwd("C:/Users/WSY/Desktop/My_data/GSEA/KEGG_Pathway_all_genes/")

FDR_genes_list = list(FDR_HS_file, FDR_NS_file, FDR_Wt_file, FDR_KO_file)
setwd("C:/Users/WSY/Desktop/My_data/GSEA/KEGG_Pathway_FDR_genes/")

SIG_genes_list = list(SIG_HS_file, SIG_NS_file, SIG_Wt_file, SIG_KO_file)
setwd("C:/Users/WSY/Desktop/My_data/GSEA/KEGG_Pathway_sig_genes/")


for (df in FDR_genes_list) {
  
  if (colnames(df)[1] == "HS_gene_id") {
    path = "./HS/"
    name = "HS"
    title_name = "KOHS vs WtHS Enriched Pathways KEGG"
    KEGG_object = "./HS/HS_KEGG.rds"
    KEGG_tsvname = "./HS/HS_KEGG.tsv"
    plot_size = c(7,5,6,4,10,5,6,6,8,8,8,6,10,6,6,6)
  } else if (colnames(df)[1] == "NS_gene_id") {
    path = "./NS/"
    name = "NS"
    title_name = "KONS vs WtNS Enriched Pathways KEGG"
    KEGG_tsvname = "./NS/NS_KEGG.tsv"
    KEGG_object = "./NS/NS_KEGG.rds"
    plot_size = c(7,5,6,4,10,5,6,6,8,8,8,6,10,6,6,6)
  } else if (colnames(df)[1] == "Wt_gene_id") {
    path = "./Wt/"
    name = "Wt"
    title_name = "WtHS vs WtNS Enriched Pathways KEGG"
    KEGG_object = "./Wt/Wt_KEGG.rds"
    KEGG_tsvname = "./Wt/Wt_KEGG.tsv"
    plot_size = c(7,5,5,4,10,5,6,6,8,6,8,6,10,6,6,6)
  } else {
    path = "./KO/"
    name = "KO"
    title_name = "KOHS vs KONS Enriched Pathways KEGG"
    KEGG_object = "./KO/KO_KEGG.rds"
    KEGG_tsvname = "./KO/KO_KEGG.tsv"
    plot_size = c(7,5,6,4,10,5,6,6,8,6,8,6,10,6,6,6)
  }
  
  print(path)
  print(name)
  print(title_name)
  
  # prepare data
  original_gene_list = df$log2FoldChange
  class(original_gene_list)
  names(original_gene_list) = df[, c(1)]
  class(names(original_gene_list))
  # View(original_gene_list)
  gene_list = na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = T)
  
  # convert gene IDs for gseKEGG function
  ids = bitr(names(original_gene_list), fromType = "ENSEMBL", 
             toType = "ENTREZID", OrgDb=organism)
  
  # remove duplicate IDs
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  remove(df2)
  df2 = df[df[, c(1)] %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # Create a vector of the gene universe
  kegg_gene_list = df2$log2FoldChange
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) = df2$Y
  
  # omit any NA values 
  kegg_gene_list = na.omit(kegg_gene_list)
  
  # sort the list in decreasing order
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  # create gseKEGG object
  kegg_organism = "mmu"
  
  kk2 = gseKEGG(geneList     = kegg_gene_list,
                organism     = kegg_organism,
                nPerm        = 10000,
                minGSSize    = 3,
                maxGSSize    = 800,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                keyType       = "ncbi-geneid")
  
  # save kk2 as csv file
  write.table(kk2, file=KEGG_tsvname, sep="\t",quote=F,row.names = F)
  KEGG_file = read.table(file=KEGG_tsvname, sep="\t", header = T)
  
  # add -log(pvalue) as a column
  kk2@result$"-log(p-value)" = -log(kk2@result$pvalue)
  kk2@result = kk2@result[order(-kk2@result$`-log(p-value)`),]
  
  # save kk2 object
  saveRDS(kk2, file = KEGG_object)
  
  # dotplot
  require(DOSE)
  dotplot_1 = dotplot(kk2, showCategory=10, title = title_name, 
                      split=".sign",label_format=30, font.size = 10) + facet_grid(.~.sign)
  pdf(paste(path, "dotplot_1_",name, ".pdf", sep = ""), width = plot_size[1], height = plot_size[2])
  print(dotplot_1)
  dev.off()
  
  dotplot_2 = dotplot(kk2, showCategory=10, title = title_name, label_format=30, font.size = 10) 
  pdf(paste(path, "dotplot_2_",name, ".pdf", sep = ""), width = plot_size[3], height = plot_size[4])
  print(dotplot_2)
  dev.off()

  # barplot
  ggdata = data.frame(value = kk2@result$`-log(p-value)`[1:10],
                      pathway = kk2@result$Description[1:10],
                      id = kk2@result$ID[1:10])
  
  barplot = ggplot(ggdata, aes(x = reorder(pathway, +value), y = value)) +    
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    labs(x = "Pathway", y = "-log(p-value)", title = title_name) +
    theme(axis.text.x = element_text(angle = 90, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10)) +
    geom_text(
      aes(label = paste(round(value, 3))),
      color = "black",size = 4,hjust=-0.1,
      position = position_dodge(0.9)) + coord_flip()
  
  pdf(paste(path, "barplot_", name, ".pdf", sep = ""), width = plot_size[5], height = plot_size[6])
  print(barplot)
  dev.off()
  
  # enrichment map
  kk3 = pairwise_termsim(kk2)
  enrichment_plot = emapplot(kk3, showCategory = 10)
  pdf(paste(path, "emapplot_", name, ".pdf", sep = ""), width = plot_size[7], height = plot_size[8])
  print(enrichment_plot)
  dev.off()
  
  # category netplot
  netplot = cnetplot(kk2, categorySize="pvalue", foldChange=gene_list, 
           showCategory = 3, cex_label_gene = 0.5)
  pdf(paste(path, "netplot_", name, ".pdf", sep = ""), width = plot_size[9], height = plot_size[10])
  print(netplot)
  dev.off()
  
  # circular centplot
  circular_cen = cnetplot(kk2, categorySize="pvalue", foldChange=gene_list, 
                          showCategory = 3, cex_label_gene = 0.5, circular = T, colorEdge = TRUE)
  pdf(paste(path, "circular_cenplot_",name, ".pdf", sep = ""), width = plot_size[11], height = plot_size[12])
  print(circular_cen)
  dev.off()

  # Ridgeplot
  library(ggridges)
  ridge = ridgeplot(kk2, showCategory = 15,
                    fill = "p.adjust",label_format = 55) + 
            labs(x = "enrichment distribution")
  pdf(paste(path, "ridgeplot_", name, ".pdf", sep = ""), width = plot_size[13], height = plot_size[14])
  print(ridge)
  dev.off()
  
  #GSEA Plot
  gsea = gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
  pdf(paste(path, "gseaplot_", name, ".pdf", sep = ""), width = plot_size[15], height = plot_size[16])
  print(gsea)
  dev.off()
  
  id_list = c(KEGG_file$ID)
  
  if (length(id_list) != 0) {
    for (id in id_list) {
      print(id)
      mmu = pathview(gene.data=kegg_gene_list, pathway.id=id, 
                     species = kegg_organism, kegg.native = T, out.suffix = name)
    }
  }
}


# Pathview ####
GO_sig_HS_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/DESeq_res_table_IPA/sig_KOHSvsWtHS_nondiscovery_10q.csv", header = T)
colnames(GO_sig_HS_file)[1] = "HS_gene_id"
GO_sig_NS_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/DESeq_res_table_IPA/sig_KONSvsWtNS_nondiscovery_10q.csv", header = T)
colnames(GO_sig_NS_file)[1] = "NS_gene_id"
GO_sig_Wt_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/DESeq_res_table_IPA/sig_WtHSvsWtNS_nondiscovery_10q.csv", header = T)
colnames(GO_sig_Wt_file)[1] = "Wt_gene_id"
GO_sig_KO_file = read.csv("C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/DESeq_res_table_IPA/sig_KOHSvsKONS_nondiscovery_10q.csv", header = T)
colnames(GO_sig_KO_file)[1] = "KO_gene_id"

sig_file_list = list(GO_sig_HS_file, GO_sig_NS_file, GO_sig_Wt_file, GO_sig_KO_file)

for (df in sig_file_list) {
  if (colnames(df)[1] == "HS_gene_id") {
    path = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/HS"
    name = "HS"
    title_name = "KOHS vs WtHS Enriched Pathways GO"
    KEGG_filename = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/HS/HS_KEGG.tsv"
    
  } else if (colnames(df)[1] == "NS_gene_id") {
    path = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/NS"
    name = "NS"
    title_name = "KONS vs WtNS Enriched Pathways GO"
    KEGG_filename = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/NS/NS_KEGG.tsv"
    
  } else if (colnames(df)[1] == "Wt_gene_id") {
    path = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/Wt"
    name = "Wt"
    title_name = "WtHS vs WtNS Enriched Pathways GO"
    KEGG_filename = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/Wt/Wt_KEGG.tsv"
    
  } else {
    path = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/KO"
    name = "KO"
    title_name = "KOHS vs KONS Enriched Pathways GO"
    KEGG_filename = "C:/Users/WSY/Desktop/My_data/IPA/no_TMM_onlyDESeq_with_raw_nondiscovery_10q/sig_GO/KO/KO_KEGG.tsv"
    
  }

  kegg_organism = "mmu"
  library(pathview)
  setwd(path)
  print(path)
  original_gene_list = df$log2FoldChange
  print(class(original_gene_list))
  names(original_gene_list) = df[, c(1)]
  print(class(names(original_gene_list)))
  
  # convert gene IDs for gseKEGG function
  ids = bitr(names(original_gene_list), fromType = "ENSEMBL", 
             toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  print("finish ids step")
  
  # remove duplication IDS
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new data frame df2 which has only the genes which were successfully mapped using the bitr function above
  remove(df2)
  df2 = df[df[, c(1)] %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # Create a vector of the gene universe
  kegg_gene_list = df2$log2FoldChange
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) = df2$Y
  
  # omit any NA values 
  kegg_gene_list = na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  # create KEGG
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")
  
  write.table(kk2,file=KEGG_filename,sep="\t",quote=F,row.names = F)
  
  KEGG_file = read.table(file=KEGG_filename, sep="\t", header = T)
  
  id_list = c(KEGG_file$ID)

  if (length(id_list) != 0) {
    for (id in id_list) {
      print(id)
      mmu = pathview(gene.data=kegg_gene_list, pathway.id=id, 
                     species = kegg_organism, kegg.native = T, out.suffix = name)
    }
  }
}