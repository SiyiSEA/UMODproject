# Part of this code used from https://github.com/meichendong/SCDC

# pre-precess scRNAseq ####
# load data
library("data.table")
library(tibble)
library(SCDC)
setwd("C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA")
file_path = "GSE107585_Mouse_kidney_single_cell_matrix/scRNA_mus_kidney_samples_1and4.csv"

# read in data
mus_count = read.csv(file_path)

# due to the file is too large, check the class of data
# if the data is tibble, than should be converted to dataframe
mus_count = as.data.frame(mus_count)
print(colnames(mus_count)[1:5])
mus_count[1:5,1:5]

# if there is a "...1" or "X" column in the header than remove it
mus_count$...1 = NULL
mus_count$X = NULL

# convert the gene name into gene id and set gene id as rowname
library("AnnotationDbi")
library(biomaRt)
library('org.Mm.eg.db')
ensembl = useMart('ensembl')
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
attributes = listAttributes(ensembl)
filters = listFilters(ensembl)
gene_name = mus_count$gene_name
info =  getBM(attributes = c("ensembl_gene_id","external_gene_name"),
              filters = "external_gene_name",
              values = gene_name,
              mart = ensembl)

# merge table based on the rows with gene name are successfully convert into gene id
merge = merge(mus_count, info, by.x = "gene_name", by.y = "external_gene_name", all.x = F, all.y = T)
merge[1:5,1:5]
merge$gene_name = NULL

# because some row with gene id are duplicated we need to remove
merge = merge[order(merge[,'ensembl_gene_id'],-merge[,'AAACCTGAGATATGCA.1']),]
merge = merge[!duplicated(merge$ensembl_gene_id),]
mus_count = merge
rownames(mus_count) = mus_count$ensembl_gene_id
mus_count$ensembl_gene_id = NULL
mus_count[1:5,1:5]
fdata = rownames(mus_count)

# read in cluster
mus_cluster = read.csv("GSE107585_Mouse_kidney_single_cell_matrix/cluster_sample_5.csv")
head(mus_cluster)
mus_cluster$...1 = NULL
mus_cluster$X = NULL
head(mus_cluster)
mus_cluster = as.data.frame(mus_cluster)

# create an SCDC rds object
mus_obj = getESET(mus_count, fdata = fdata, pdata = mus_cluster)
rownames(mus_obj)[1:4]
saveRDS(mus_obj, file = "scRNA_mus_kidney_5_geneid.rds")

# demoplot
color_1 = c("#f39595", "#ea437d","#d56cef", "#cdbff8",
          "#6c9eef", "#43d4ea", "#95f3d7", "#43ea61", 
          "#109328", "#c0e619","#bc7915","#faa307",
          "#c4c4c4", "#653e3e", "#d28989", "#3b3b3b")

interest_name = c("Endo", "Podo", "PT", "LOH", 
                  "DCT", "CD-PC", "CD-IC", "CD-Tran", 
                  "Novel1", "Fib", "Marco", "Neutro", 
                  "B-lymph", "T-lymph","NK","Novel2")

pdf(file = "C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/plot1_SCDC/Demoplot_Sample_5.pdf", width = 10, height = 8)
demo = DemoPlot(mus_obj, cluster = "Cluster_Name", sample = "Sample", select.ct = interest_name, Palette = color_1)
dev.off()

# qc for one subject in the object
mus_obj.qc = SCDC_qc_ONE(mus_obj, ct.varname = "Cluster_Name", sample = "Sample", scsetname = "Mus_kidney",
                         ct.sub = interest_name, qcthreshold = 0.7, cbPalette = color) 

# qc for more than one subjects in the object
mus_obj.qc = SCDC_qc(mus_obj, ct.varname = "Cluster_Name", sample = "Genotype", scsetname = "Mus_kidney",
                         ct.sub = interest_name, qcthreshold = 0.7, cbPalette = color)

saveRDS(mus_obj.qc, file = "scRNA_mus_kidney_5_qc_geneid.rds")

# heatmap 
file_name = "C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/plot1_SCDC/heat_Samples_5.pdf"
pdf(file = file_name, width = 10, height = 8)
heatplot = mus_obj.qc$heatfig
print(heatplot)
dev.off()

# pre-recess bulk RNA ####
nondis = read.csv("gene_count_matrix_nondiscovery_10q.csv")
row.names(nondis) = nondis$gene_id
nondis$gene_id = NULL
fdata_nondis = rownames(nondis)
pdata_nondis = read.csv("experiment_design_table.csv")
nondiseset = getESET(nondis, fdata = fdata_nondis, pdata = pdata_nondis)
saveRDS(nondiseset, file = "bulkRNA_UMOD.rds")
class(eset)


# SCDC for one subject ####
# read in data ####
library(SCDC)
sample_5_qc <- readRDS("C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/scRNA_mus_kidney_5_qc_geneid.rds")
sample_6_qc <- readRDS("C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/scRNA_mus_kidney_6_qc_geneid.rds")

sample_5 = sample_6_qc$sc.eset.qc
sample_6 = sample_6_qc$sc.eset.qc

bulkRNA_UMOD <- readRDS("C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/bulkRNA_UMOD.rds")

# do ensemble on each groups
# choose sample 5 and sample 6 reference scRNA
sc_ref = list(ref5 = sample_5_qc$sc.eset.qc, ref6 = sample_6_qc$sc.eset.qc)
all_clusters = c("Endo", "Podo", "PT", "LOH", 
                 "DCT", "CD-PC", "CD-IC", "CD-Tran", 
                 "Novel1", "Fib", "Marco", "Neutro", 
                 "B-lymph", "T-lymph","NK","Novel2")

bulkRNA_UMOD_KOHS_ens = SCDC_ENSEMBLE(bulk.eset = bulkRNA_UMOD[, bulkRNA_UMOD$Condition == "KOHS"],
                                      sc.eset.list = sc_ref,
                                      ct.varname = "Cluster_Name", sample = "Sample", 
                                      ct.sub = all_clusters, grid.search = T, search.length = 0.5)

bulkRNA_UMOD_KONS_ens = SCDC_ENSEMBLE(bulk.eset = bulkRNA_UMOD[, bulkRNA_UMOD$Condition == "KONS"],
                                      sc.eset.list = sc_ref, 
                                      ct.varname = "Cluster_Name", sample = "Sample", 
                                      ct.sub = all_clusters, grid.search = T, search.length = 0.5)

bulkRNA_UMOD_WtHS_ens = SCDC_ENSEMBLE(bulk.eset = bulkRNA_UMOD[, bulkRNA_UMOD$Condition == "WtHS"],
                                      sc.eset.list = sc_ref, 
                                      ct.varname = "Cluster_Name", sample = "Sample", 
                                      ct.sub = all_clusters, grid.search = T, search.length = 0.5)

bulkRNA_UMOD_WtNS_ens = SCDC_ENSEMBLE(bulk.eset = bulkRNA_UMOD[, bulkRNA_UMOD$Condition == "WtNS"],
                                      sc.eset.list = sc_ref, 
                                      ct.varname = "Cluster_Name", sample = "Sample", 
                                      ct.sub = all_clusters, grid.search = T, search.length = 0.5)

# plot proprotion boxplot on (KOHS vs WtHS) samples ####
library(reshape2)
library(ggplot2)
# getPropBox function
getPropBox <- function(ens_h, ens_d, metric = NULL, ref, input_wt = NULL){
  if (!is.null(input_wt)){
    prop.h = as.data.frame(wt_prop(input_wt, ens_h$prop.only))
    prop.d = as.data.frame(wt_prop(input_wt, ens_d$prop.only))
  } else {
    prop.h = as.data.frame(wt_prop(ens_h$w_table[metric,1:2], ens_h$prop.only))
    prop.d = as.data.frame(wt_prop(ens_d$w_table[metric,1:2], ens_d$prop.only))
  }
  prop2 <- rbind(prop.h, prop.d)
  prop2$condition <- c(rep("WtHS", nrow(prop.h)), rep("WtNS", nrow(prop.d)))
  # HS group: "KOHS" "WtHS"
  # NS group: "KONS" "WtNS"
  # KO group: "KOHS" "KONS"
  # Wt group: "WtHS" "WtNS"
  dtmelt <- reshape2::melt(prop2, id.vars = "condition")
  dtmelt$ref <- as.factor(ref)
  return(dtmelt)
}

# prepare data for PropBox
# HS group: "KOHS" "WtHS" ens_h = bulkRNA_UMOD_KOHS_ens, ens_d = bulkRNA_UMOD_WtHS_ens,
# NS group: "KONS" "WtNS" ens_h = bulkRNA_UMOD_KONS_ens, ens_d = bulkRNA_UMOD_WtNS_ens,
# KO group: "KOHS" "KONS" ens_h = bulkRNA_UMOD_KOHS_ens, ens_d = bulkRNA_UMOD_WtHS_ens,
# Wt group: "WtHS" "WtNS" ens_h = bulkRNA_UMOD_WtHS_ens, ens_d = bulkRNA_UMOD_WtNS_ens,

ens_spearman <- getPropBox(ens_h = bulkRNA_UMOD_WtHS_ens, ens_d = bulkRNA_UMOD_WtNS_ens,
                              metric = 1, 
                              ref = "ENSEMBLE+SpearmanR")
ens_sample1 <- getPropBox(ens_h = bulkRNA_UMOD_WtHS_ens, ens_d = bulkRNA_UMOD_WtNS_ens,
                             ref = "scRNA-sample-5", 
                             input_wt = c(1,0))
ens_sample2 <- getPropBox(ens_h = bulkRNA_UMOD_WtHS_ens, ens_d = bulkRNA_UMOD_WtNS_ens,
                             ref = "scRNA-sample-6", 
                             input_wt = c(0,1))

dtall <- rbind(ens_spearman, ens_sample1, ens_sample2)
dtall$refcond <- paste(dtall$ref, dtall$condition)

colfunc <- colorRampPalette(c("red", "white"))
colfunc2 <- colorRampPalette(c("blue", "white"))

# plot ProBox
# ref12 HS group
#ProBox_label = c("ENSEMBLE+SpearmanR KOHS", "scRNA-sample-1 KOHS",
#                 "scRNA-sample-2 KOHS", "ENSEMBLE+SpearmanR WtHS",
#                 "scRNA-sample-1 WtHS","scRNA-sample-2 WtHS")

# HS group
ProBox_label = c("ENSEMBLE+SpearmanR KOHS", "scRNA-sample-5 KOHS",
                 "scRNA-sample-6 KOHS", "ENSEMBLE+SpearmanR WtHS",
                 "scRNA-sample-5 WtHS","scRNA-sample-6 WtHS")
#NS group
ProBox_label = c("ENSEMBLE+SpearmanR KONS", "scRNA-sample-5 KONS",
                 "scRNA-sample-6 KONS", "ENSEMBLE+SpearmanR WtNS",
                 "scRNA-sample-5 WtNS","scRNA-sample-6 WtNS")
# KO group
ProBox_label = c("ENSEMBLE+SpearmanR KOHS", "scRNA-sample-5 KOHS",
                 "scRNA-sample-6 KOHS", "ENSEMBLE+SpearmanR KONS",
                 "scRNA-sample-5 KONS","scRNA-sample-6 KONS")

# Wt group
ProBox_label = c("ENSEMBLE+SpearmanR WtHS", "scRNA-sample-5 WtHS",
                 "scRNA-sample-6 WtHS", "ENSEMBLE+SpearmanR WtNS",
                 "scRNA-sample-5 WtNS","scRNA-sample-6 WtNS")


pfa2 <- ggplot(dtall[dtall$refcond %in% ProBox_label,], 
               aes(x=variable, y=value, 
                   color = factor(refcond, levels=ProBox_label))) + 
  geom_boxplot(outlier.size=-1) +
  geom_jitter(aes(x=variable,
                  color = factor(refcond, levels=ProBox_label)),
              position = position_dodge(0.75), alpha=0.5,cex=0.25)+
  theme(axis.text.x = element_text(angle = 30, hjust=1, size=10),
        axis.text.y = element_text(size = 10),
        text = element_text(size = 10),
        plot.title = element_text(size=10, face = "bold"),
        plot.margin=unit(c(1,1,-5,0), "mm"),
        legend.position="top",legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box.spacing = unit(0, "mm"))+
  scale_color_manual(values=c(colfunc(10)[seq(1,9,3)],colfunc2(10)[seq(1,9,3)])) +
  ylab("") + xlab("")
pfa2

# save proBox plot
pdf(file = "C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/plot1_SCDC/PropBox_ref56_Wt.pdf",width = 10, height = 6)
print(pfa2)
dev.off()


# visualizes the SCDC estimation results ####
# prepare data
# take KOHS as an example
# KOHS bulk obj = bulkRNA_UMOD_KOHS_ens
library(ggplot2)
library(reshape2)

combo_WtHS <- lapply(bulkRNA_UMOD_WtHS_ens$prop.list, function(x){
  x$prop.est[,all_clusters]
})
combo_WtHS$ref5 = colMeans(combo_WtHS$ref5)
combo_WtHS$ref6 = colMeans(combo_WtHS$ref6)

#stimulate scRNA proportation by Seurat
# we can use one of the scRNA ref sample
sample_5_pseudo_all <- generateBulk_allcells(sample_5, 
                                             ct.varname = "Cluster_Name", 
                                             sample = 'Sample', 
                                             ct.sub = all_clusters, 
                                             disease = "Genotype")

# getPearson: calculating the Pearson correlation
getPearson <- function(dt1, dt2){
  cor(dt1, dt2)
}
# ens_res:summarize the ENSEMBLE results
ens_res <- function(wt, proplist, pseudo_all){
  ens <- wt_prop(wt, proplist = proplist) 
  p.ens <- getPearson(pseudo_all$truep, ens)
  return(list(prop = ens, p = p.ens))
}

getENSplot <- function(enstable, combo, subject, subcat = NULL, pseudo_all, file_name){
  if (!is.null(subcat)){
    enstable <- enstable[subcat,]
  }
  wtres.list <- list()
  
  for (i in 1:nrow(enstable)){
    wtres.list[[i]] <- ens_res(wt = enstable[i,1:2], proplist = combo, pseudo_all = pseudo_all)$prop
  }
  
  names(wtres.list) <- rownames(enstable)
  pcor.list <- NULL
  for (i in 1:nrow(enstable)){
    pcor.list[[i]] <- ens_res(wt = enstable[i,1:2], proplist = combo, pseudo_all = pseudo_all)$p
  }
  
  names(pcor.list) <- rownames(enstable)
  est.all <- rbind(pseudo_all$truep,
                   do.call(rbind, combo), 
                   do.call(rbind, wtres.list))
  
  rownames(est.all) <- as.factor(c("Seurat", "scRNA ref1", "scRNA ref2 R", 
                                   paste("ENSEMBLE:",names(pcor.list))))
  
  p_t <- getPearson(pseudo_all$truep, combo[[1]])
  p_p <- getPearson(pseudo_all$truep, combo[[2]])
  dat_text <- data.frame(label = c("",paste("R=", round(c(p_t, p_p, unlist(pcor.list)),2), sep = "")),
                         Var1 = as.factor(c("Seurat", "scRNA ref1", "scRNA ref2", 
                                            paste("ENSEMBLE:",names(pcor.list)))),
                         value =rep(0.9, length(pcor.list)+3),
                         Var2 = rep("Podo",length(pcor.list)+3))
  meltdt <- reshape2::melt(est.all)
  P <- ggplot(meltdt, aes(x=Var1, y=value, fill = factor(Var2, levels = all_clusters))) +
    geom_bar(stat = "identity") + xlab("") + ylab("Percentage") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size=10),
          axis.text.y = element_text(size = 10),
          text = element_text(size = 10),
          plot.title = element_text(size=10, face = "bold"),
          plot.margin=unit(c(1,1,-5,0), "mm"),
          legend.position="top",legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.box.spacing = unit(0, "mm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    geom_text(data = dat_text,
              mapping = aes(x=Var1, y=value,label = label),
              hjust   = 0.5, vjust   = -3.5) + 
    scale_fill_manual(values=c("#f39595", "#ea437d","#d56cef", "#cdbff8",
                               "#6c9eef", "#43d4ea", "#95f3d7", "#43ea61", 
                               "#109328", "#c0e619","#bc7915","#faa307",
                               "#c4c4c4", "#653e3e", "#d28989", "#3b3b3b"))
  P
  pdf(file = "C:/Users/WSY/Desktop/My_data/SCDC/mouse_kidney_scRNA/plot1_SCDC/decon_ref56_KOHS_difcor.pdf", width = 6, height = 5)
  print(P)
  dev.off()
}


getENSplot(enstable = bulkRNA_UMOD_WtHS_ens$w_table, 
           combo=combo_WtNS, 
           subcat = c("Spearman_Y","LAD"), 
           pseudo_all=sample_5_pseudo_all)
