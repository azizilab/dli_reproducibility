rm(list=ls())

### set directory
setwd("/Users/jamesbrenner/My Drive/final_data/km_update/")

#install.packages("remotes")
#remotes::install_github("satijalab/seurat-disk")
#remotes::install_github("mojaveazure/seurat-disk")

library(dplyr)
library(Seurat)
#library(SeuratData)
#library(SeuratDisk)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(gapminder)
library(stringr)
library(Matrix)
dir.create("plots")
dir.create("objects")

### Importing data
cell_matrix <- fread(input = "aml_X.csv", sep = ',', header = TRUE)
gene_names <- read.csv("aml_obs.csv", header = FALSE)
cell_bar <- read.csv("aml_var.csv", header = FALSE)
cell_matrix <- cell_matrix[,-1]
#colnames(cell_matrix) <- as.character(cell_bar)
#rownames(cell_matrix) <- gene_names
colnames(cell_matrix) <- cell_bar$V1
rownames(cell_matrix) <- gene_names$V1

seurat_object <- CreateSeuratObject(counts = cell_matrix)


all_cells <- ReadMtx(mtx = "aml_X.csv", 
        cells = "aml_var.csv",     
        features = "aml_obs.csv",  
        cell.column = 1,
        feature.column = 1,      
        cell.sep = ",",
        feature.sep = ",",
        skip.cell = 0,
        skip.feature = 0,
        mtx.transpose = TRUE,    
        unique.features = TRUE,
        strip.suffix = FALSE)


seurat_object <- CreateSeuratObject(counts = cell_matrix, row.names = gene_names)
test_cells <- seurat_object
meta <- data.frame( fread('allBM_metadata.csv') )
rownames(meta) <- meta$barcode
test_cells@meta.data <- meta

#test_cells$random <- Idents(test_cells)
### Subset T and NK cells 
Idents(test_cells) <- "cluster_number"
tcell.clusters <- c("0", "3", "7", "17", "31", "37", "45")
tcell.clusters.cd8 <- c("0", "7", "31", "37", "45")
nk.clusters <- c("2", "15", "30")
pbmc_tcell_nk <- subset(test_cells, idents = c(tcell.clusters, nk.clusters))

### Scaling data
#data is normalized, log transformed, but not scaled 
pbmc_tcell_nk <- FindVariableFeatures(pbmc_tcell_nk, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell_nk)
pbmc_tcell_nk <- ScaleData(pbmc_tcell_nk) 
pbmc_tcell_nk <- RunPCA(pbmc_tcell_nk, features = VariableFeatures(object = pbmc_tcell_nk))
pbmc_tcell_nk <- FindNeighbors(pbmc_tcell_nk, dims = 1:10)
pbmc_tcell_nk <- FindClusters(pbmc_tcell_nk, resolution = 1)

# Run UMAP/tSNE
pbmc_tcell_nk <- RunUMAP(pbmc_tcell_nk, dims = 1:10)
Idents(pbmc_tcell_nk) <- "cluster_number"
DimPlot(pbmc_tcell_nk, reduction = "umap", label = TRUE)
FeaturePlot(pbmc_tcell_nk, features = c("ZNF683", "CD8A", "CD4", "IL2RB"))

saveRDS(pbmc_tcell_nk, file = "objects/km_tcell_nk.rds")

gene.list.nk <- c("ZNF683", "CD8A", "CD4", "CD3E", "BCL11B", 
                  "SYK", "KIR2DL3", "NCR1", 
                  "GZMB", "KIR2DL1", "KIR2DL4", "KIR3DL2",
                  "KLRC2", "KLRC3", "TYROBP", "CD244", 
                  "KLRC1", "KIR3DL1", "PLCG2", "FCGR3A",
                  "ITGAX", "S1PR5", "PRF1", "TBX21", "GNLY",
                  "IL2RB", "LAT2", "LYN", "NCAM1", "SLAMF7", 
                  "TCF7", "CD27", "CD5", "PDCD1", "CD28", 
                  "CCR7", "CD6", "IL7R", "THEMIS", "TESPA1", 
                  "B3GAT1", "XCL1", "XCL2", "CD160", "SELL", "CXCR6")
# Following lines are used to format for plotting
num_idents_test <- length(unique(Idents(object = pbmc_tcell_nk)))
na_colors_test <- rep("NA", num_idents_test)
pbmc_tcell_nk_test <- pbmc_tcell_nk
my_levels <- c("2",  "15", "30", "0", "3", "7", "17", "31", "37", "45")
# Re-level object@ident
pbmc_tcell_nk_test@active.ident <- factor(x = pbmc_tcell_nk_test@active.ident, levels = my_levels)
combined_averages <- AverageExpression(pbmc_tcell_nk_test, return.seurat = TRUE) 
r_avg <- AverageExpression(subset(pbmc_tcell_nk_test, subset = response == "RESPONDER"), return.seurat = TRUE)
nr_avg <- AverageExpression(subset(pbmc_tcell_nk_test, subset = response == "NONRESPONDER"), return.seurat = TRUE)

### Plot heat maps 
nk <- DoHeatmap(combined_averages, features = gene.list.nk, group.by = "ident", label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "All") 
nk
ggsave(file = "plots/km_all_tcell_heat_nk_panel.eps", width = 7, height = 8, dpi = 150,)

r.heat <- DoHeatmap(r_avg, features = gene.list.nk, group.by = "ident", label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "R") 
r.heat
nr.heat <- DoHeatmap(nr_avg, features = gene.list.nk, group.by = "ident", label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "NR")
r.heat + nr.heat
ggsave(file = "plots/km_all_tcell_heat_nk_panel_r_nr.eps", width = 11, height = 8, dpi = 150)

### Subset only T cells (CD4 and CD8) and scale data
pbmc_tcell_only <- subset(test_cells, idents = tcell.clusters)
pbmc_tcell_only <- FindVariableFeatures(pbmc_tcell_only, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell_only)
pbmc_tcell_only <- ScaleData(pbmc_tcell_only) 
pbmc_tcell_only <- RunPCA(pbmc_tcell_only, features = VariableFeatures(object = pbmc_tcell_only))
pbmc_tcell_only <- FindNeighbors(pbmc_tcell_only, dims = 1:10)
pbmc_tcell_only <- FindClusters(pbmc_tcell_only, resolution = 1)

# Run UMAP/tSNE
pbmc_tcell_only <- RunUMAP(pbmc_tcell_only, dims = 1:10)
Idents(pbmc_tcell_only) <- "cluster_number"
DimPlot(pbmc_tcell_only, reduction = "umap", label = TRUE)
FeaturePlot(pbmc_tcell_only, features = c("ZNF683", "CD8A", "CD4", "IL2RB"))

saveRDS(pbmc_tcell_only, file = "objects/km_tcell.rds")

my_levels_2 <- c("0", "3", "7", "17", "31", "37", "45") 
pbmc_tcell_only_test <- pbmc_tcell_only

# Re-level object@ident
Idents(pbmc_tcell_only_test) <- "cluster_number"
pbmc_tcell_only_test@active.ident <- factor(x = pbmc_tcell_only_test@active.ident, levels = my_levels_2)
combined_averages_tcell <- AverageExpression(pbmc_tcell_only_test, return.seurat = TRUE) 


gene.list.tcell <- c("SYNE1", "ICAM1", "CCL4L2", "GNLY", "DTHD1", 
                 "NKG7", "PLAUR", "IL1RN", "ANPEP", "CCL2", "MMP19",
                 "SERPINB2", "SELL", "DDX3Y", "DLG5", "SETBP1", "GZMH", 
                 "PZP", "CXCL3", "IL1B", "LRRN3",
                 "CCL4", "GZMB", "XCL2", "PRF1", "IFNG", "CXCR6", "KIR3DL1")

### Dot and heat plots of features/panels 
t <- DoHeatmap(combined_averages_tcell, features = gene.list.tcell, group.by = "ident", label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "All") 
t

gene.list.tcell.cd8 <- c("CD8A", "CD4", "CD69", "ITGAE", "KLRC2", "CD6", "ENTPD1", 
                         "TCF7", "LEF1", "CCR7", "SELL", "MAL", 
                         "IL7R", "GPR183", "ZFP36L2", "CXCR4", 
                         "ZNF683", "CD52", "HOPX", "ID2", "CXCR6", "XCL1", "XCL2", 
                         "TBX21", "ASCL2", "CX3CR1", "KLRG1", 
                         "KLRD1", "TYROBP", "KIR2DL3", "KIR2DL1", "KIR3DL1", "KIR3DL2", "CD160", "EOMES", "TXK", "KLRC1", "KIR2DL4",
                         "GZMK", "CXCR5", "CCR4", "CD28", "CXCR3", "GZMH", "CD27", "HLA-DRB1", 
                         "PDCD1", "CXCL13", "LAYN", 
                         "STAT1", "IFIT1", "ISG15", "CCR1", 
                         "SLC4A10", "KLRB1", "TMIGD2", "RORA", "RORC", "ZBTB16", "IL26", "IL23R", 
                         "NME1", "NME2", "MND1", "SPC24", "MYB")

gene.list.tcell.cd4 <- c("CD40LG", 
                         "FOXP3", 
                         "TCF7", "LEF1", "TXK", "CCR7", "SELL", "MAL", "CXCR5", "ADSL", "IL16", "IL7R", 
                         "TNF", "AREG", "TIMP1", "CREM", "CCL5", "CAPG", "GZMK", "KLRG1", "CX3CR1", "TBX21", 
                         "RORA", "RORC", "CCR6", "IL23R", "IL26", 
                         "TOX", "TOX2", "IL21", "CXCL13", "GNG4", "CD200", "BCL6", "ZBED2", "CCL3", "CCL4", "IFNG", "GZMB", "LAG3", "HAVCR2", 
                         "RTKN2", "IL2RA", "S1PR1", "TNFRSF9", "CTLA4", "LAYN", "STAT1", "IFIT1", "IRF7", 
                         "NME1", "NME2", "MND1", "SPC24", "CCR4")
# This is from giacomo's paper 
gene.list.tcell.num1 <- c("CD3E", "CD4", "CD8A", 
                          "SELL", "CCR7", "IL7R", "CD28", "FAS", "CD27", "ITGAE", "ITGAL", "ITGAM", "ITGAX", 
                          "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "CD244", "KLRG1", "TNFRSF14", "BTLA", "CD160", 
                          "CD38", "ENTPD1", "NT5E", "CD69", "IL2RA", "ICOS", "TNFRSF4", "TNFRSF9", "HLA-DRA", "CD40LG", 
                          "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "IFNG", "FASLG", "TNF", "IL2",  
                          "LEF1", "TCF7", "EOMES", "TBX21", "PRDM1", "TOX", "GATA3", "ID2", "ID3", "NR4A1", "ZNF683", "FOXP3", "MKI67", "TOP2A",  
                          "TRGV9", "TRDV2", "KLRB1", "KLRC3")

gene.list.tcell.act.exh <- c("ENTPD1", "IL2RA", "ICOS", "TNFRSF4", "CD40LG",
                             "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4")

DotPlot(object = pbmc_tcell_only_test, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_all_tcell_dot_cd8.eps", width = 15, height = 8, dpi = 150,)
    
DotPlot(object = pbmc_tcell_only_test, features = gene.list.tcell.cd4, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_all_tcell_cd4.eps", width = 15, height = 8, dpi = 150,)

DotPlot(object = pbmc_tcell_only_test, features = gene.list.tcell.num1, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 

DoHeatmap(combined_averages_tcell, features = gene.list.tcell.num1, group.by = "ident", label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "Giacomo Paper") 
ggsave(file = "plots/km_all_tcell_giacomo.eps", width = 15, height = 8, dpi = 150)

### Subset only cluster 0 and re-cluster
pbmc_tcell_0 <- subset(test_cells, idents = "0")
#data is normalized, log transformed, but not scaled 
pbmc_tcell_0 <- FindVariableFeatures(pbmc_tcell_0, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell_0)
pbmc_tcell_0 <- ScaleData(pbmc_tcell_0) 
pbmc_tcell_0 <- RunPCA(pbmc_tcell_0, features = VariableFeatures(object = pbmc_tcell_0))
pbmc_tcell_0 <- FindNeighbors(pbmc_tcell_0, dims = 1:10)
pbmc_tcell_0 <- FindClusters(pbmc_tcell_0, resolution = 1)

# Run UMAP/tSNE
pbmc_tcell_0 <- RunUMAP(pbmc_tcell_0, dims = 1:10)
Idents(pbmc_tcell_0) <- "seurat_clusters"
#Idents(pbmc_tcell_0) <- "response"
#Idents(pbmc_tcell_0) <- "time"
DimPlot(pbmc_tcell_0, reduction = "umap", label = FALSE)
#ggsave(file = "plots/km_cluster0_umap_clusters.eps")

saveRDS(pbmc_tcell_0, file = "objects/km_cluster0.rds")

FeaturePlot(pbmc_tcell_0, features = c("ZNF683", "CD8A", "CD4", "CD3E", 
                                       "TIGIT", "CD226", "GZMH", "IL2RB", 
                                       "KLRG1", "CX3CR1", "IL7R", "FCGR3A", 
                                       "TBX21", "BACH2"))
ggsave(file = "plots/km_cluster0_umap_features.eps", width = 12, height = 12, dpi = 150)

### Dot plot
#pbmc_tcell_0@active.ident <- factor(x = pbmc_tcell_0@active.ident, levels = my_levels_2)
combined_averages_0 <- AverageExpression(pbmc_tcell_0, return.seurat = TRUE) 
DotPlot(object = pbmc_tcell_0, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_cluster0_dot_cd8_panel.eps", width = 15, height = 8, dpi = 150)

r_avg_0 <- AverageExpression(subset(pbmc_tcell_0, subset = response == "RESPONDER"), return.seurat = TRUE)
nr_avg_0 <- AverageExpression(subset(pbmc_tcell_0, subset = response == "NONRESPONDER"), return.seurat = TRUE)

DotPlot(object = subset(pbmc_tcell_0, subset = response == "RESPONDER"), features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) +
  labs(title = "R") 
ggsave(file = "plots/km_cluster0_dot_cd8_panel_r.eps", width = 15, height = 8, dpi = 150)

DotPlot(object = subset(pbmc_tcell_0, subset = response == "NONRESPONDER"), features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) +
  labs(title = "NR") 
ggsave(file = "plots/km_cluster0_dot_cd8_panel_nr.eps", width = 15, height = 8, dpi = 150)

### Subset only CD8 T cells 
pbmc_tcell_cd8 <- subset(test_cells, idents = tcell.clusters.cd8)
#data is normalized, log transformed, but not scaled 
pbmc_tcell_cd8 <- FindVariableFeatures(pbmc_tcell_cd8, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell_cd8)
pbmc_tcell_cd8 <- ScaleData(pbmc_tcell_cd8) 
pbmc_tcell_cd8 <- RunPCA(pbmc_tcell_cd8, features = VariableFeatures(object = pbmc_tcell_cd8))
pbmc_tcell_cd8 <- FindNeighbors(pbmc_tcell_cd8, dims = 1:10)
pbmc_tcell_cd8 <- FindClusters(pbmc_tcell_cd8, resolution = 1)

# Run UMAP/tSNE
# Creating new column with cluster and response for plotting
pbmc_tcell_cd8$ID <- paste(pbmc_tcell_cd8$response, pbmc_tcell_cd8$cluster_number, sep = "_")
#pbmc_tcell_cd8["ID"] <- str_replace(pbmc_tcell_cd8[['ID']], "NONRESPONDER", "NR")
pbmc_tcell_cd8 <- RunUMAP(pbmc_tcell_cd8, dims = 1:10)
Idents(pbmc_tcell_cd8) <- "ID" 
DimPlot(pbmc_tcell_cd8, reduction = "umap", label = FALSE)

saveRDS(pbmc_tcell_cd8, file = "objects/km_cd8s.rds")

my_levels_3 <- c("RESPONDER_0", "NONRESPONDER_0", "RESPONDER_7", "NONRESPONDER_7",              
                 "RESPONDER_31", "NONRESPONDER_31", "RESPONDER_37", "NONRESPONDER_37", 
                 "RESPONDER_45", "NONRESPONDER_45") 
Idents(pbmc_tcell_cd8) <- "ID" 
pbmc_tcell_cd8@active.ident <- factor(x = pbmc_tcell_cd8@active.ident, levels = my_levels_3)


pbmc_tcell_cd8_avg <- AverageExpression(pbmc_tcell_cd8, return.seurat = TRUE) 
DoHeatmap(pbmc_tcell_cd8_avg, features = gene.list.tcell.cd8, group.by = "ident", label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "R/NR CD8") 
ggsave(file = "plots/km_cd8only_heat_cd8_panel_r_nr.eps", width = 8, height = 12, dpi = 150)

DotPlot(object = pbmc_tcell_cd8, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=135)) +
  theme(axis.text.y = element_text(hjust=0.01)) +
  theme(axis.text.y = element_text(vjust=0.01)) +
  labs(title = "R/NR CD8") 
ggsave(file = "plots/km_cd8only_dot_cd8_panel_r_nr.eps", width = 15, height = 8, dpi = 150)

Idents(pbmc_tcell_cd8) <- "cluster_number" 
pbmc_tcell_cd8@active.ident <- factor(x = pbmc_tcell_cd8@active.ident, levels = my_levels_2)
combined_averages_cd8 <- AverageExpression(pbmc_tcell_cd8, return.seurat = TRUE) 
DotPlot(object = pbmc_tcell_cd8, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_cd8only_dot_cd8_panel.eps", width = 15, height = 8, dpi = 150)

DotPlot(object = pbmc_tcell_cd8, features = gene.list.tcell.act.exh, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_cd8only_dot_act_exh_panel.eps", width = 6, height = 4, dpi = 150)


#### Subset cells based on hobit expression 

pbmc_tcell_hobithigh <- subset(pbmc_tcell_0, subset = ZNF683 > 1)
pbmc_tcell_hobitlow <- subset(pbmc_tcell_0, subset = ZNF683 < 1)
pbmc_tcell_hobithigh <- FindNeighbors(pbmc_tcell_hobithigh, dims = 1:10)
pbmc_tcell_hobithigh <- FindClusters(pbmc_tcell_hobithigh, resolution = 1)
pbmc_tcell_hobitlow <- FindNeighbors(pbmc_tcell_hobitlow, dims = 1:10)
pbmc_tcell_hobitlow <- FindClusters(pbmc_tcell_hobitlow, resolution = 1)

pbmc_tcell_hobithigh_test <- pbmc_tcell_hobithigh
pbmc_tcell_hobithigh_test$hobit <- "high"
pbmc_tcell_hobitlow_test <- pbmc_tcell_hobitlow
pbmc_tcell_hobitlow_test$hobit <- "low"
pbmc_tcell_hobit_comb <- merge(pbmc_tcell_hobithigh_test, y = pbmc_tcell_hobitlow_test, project = "new_object")


# Run UMAP/tSNE
pbmc_tcell_hobithigh <- RunUMAP(pbmc_tcell_hobithigh, dims = 1:10)
Idents(pbmc_tcell_hobithigh) <- "seurat_clusters"
DimPlot(pbmc_tcell_hobithigh, reduction = "umap", label = FALSE)
#ggsave(file = "km_cluster0_umap_hobithigh.eps")
FeaturePlot(pbmc_tcell_hobithigh, features = c("ZNF683", "CD8A", "CD4", "CD3E", 
                                               "IL7R", "TBX21", "CX3CR1", "GZMH", "KLRG1", 
                                               "KLRD1", "KIR2DL3", "KIR2DL1", "KIR3DL1", "KIR3DL2", "KLRC1", "KLRC2"))
ggsave(file = "plots/km_cluster0_umap_features_hobithigh.eps", width = 12, height = 12, dpi = 150)

pbmc_tcell_hobitlow <- RunUMAP(pbmc_tcell_hobitlow, dims = 1:10)
Idents(pbmc_tcell_hobitlow) <- "seurat_clusters"
DimPlot(pbmc_tcell_hobitlow, reduction = "umap", label = FALSE)
#ggsave(file = "plots/km_cluster0_umap_hobitlow.eps")
FeaturePlot(pbmc_tcell_hobitlow, features = c("ZNF683", "CD8A", "CD4", "CD3E", 
                                       "IL7R", "TBX21", "CX3CR1", "GZMH", "KLRG1", 
                                       "KLRD1", "KIR2DL3", "KIR2DL1", "KIR3DL1", "KIR3DL2", "KLRC1", "KLRC2"))
ggsave(file = "plots/km_cluster0_umap_features_hobitlow.eps", width = 12, height = 12, dpi = 150)

combined_averages_hobithigh <- AverageExpression(pbmc_tcell_hobithigh, return.seurat = TRUE) 
combined_averages_hobitlow <- AverageExpression(pbmc_tcell_hobitlow, return.seurat = TRUE) 
DotPlot(object = pbmc_tcell_hobithigh, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_cluster0_hobithigh.eps", width = 15, height = 8, dpi = 150)

DotPlot(object = pbmc_tcell_hobitlow, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/km_cluster0_hobitlow.eps", width = 15, height = 8, dpi = 150)


## Looking at dot plots for hobit high/low - binary 
Idents(pbmc_tcell_hobit_comb) <- "hobit"
combined_averages_hobit <- AverageExpression(pbmc_tcell_hobit_comb, return.seurat = TRUE) 
DotPlot(object = pbmc_tcell_hobit_comb, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 


violin <- c("ZNF683", "CD8A", "CD4", "CD3E", 
            "IL7R", "TBX21", "CX3CR1", "GZMH", "KLRG1", 
            "KLRD1", "KIR2DL3", "KIR2DL1", "KIR3DL1", "KIR3DL2", "KLRC1", "KLRC2")
Idents(pbmc_tcell_hobit_comb) <- "hobit"
VlnPlot(pbmc_tcell_hobit_comb, features = gene.list.tcell.cd8, slot = "counts", log = TRUE)
ggsave(file = "plots/km_cluster0_hobit_violin.png", width = 15, height = 45, dpi = 100)

hobit.markers <- FindMarkers(pbmc_tcell_hobit_comb, ident.1 = "high", ident.2 = "low", min.pct = 0.1, logfc.threshold = 0.1)
head(hobit.markers, n = 92)
VlnPlot(pbmc_tcell_hobit_comb, features = c("ZNF683", "CMC1", "TRAV19"), slot = "counts", log = TRUE)
ggsave(file = "plots/km_cluster0_hobit_violin_diff.png")

## QC on each hobit seurat object before merging 
VlnPlot(pbmc_tcell_hobithigh, features = gene.list.tcell.cd8, slot = "counts", log = TRUE)
ggsave(file = "plots/testhigh.png", width = 15, height = 45, dpi = 100)
VlnPlot(pbmc_tcell_hobitlow, features = gene.list.tcell.cd8, slot = "counts", log = TRUE)
ggsave(file = "plots/testlow.png", width = 15, height = 45, dpi = 100)

