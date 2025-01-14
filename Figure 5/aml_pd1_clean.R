rm(list=ls())

### set directory
setwd("")

library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(gapminder)
library(ggpubr)
dir.create("plots")
dir.create("objects")

# Convert downloaded files into something readable by seurat 
load("GSM5936942_readCounts_TME_NL.rda")
meta <- data.frame(fread("GSM5936942_meta_data_TME_NL.txt.gz"))
AML <- CreateSeuratObject(counts = readCounts, project = "AML_Immune", min.cells = 3)
rownames(meta) <- meta$V1
AML@meta.data <- meta
pbmc <- AML

# Load the PBMC dataset
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Take a subset and normalize data
pbmc_sub <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc_sub <- NormalizeData(pbmc_sub, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_sub <- FindVariableFeatures(pbmc_sub, selection.method = "vst", nfeatures = 2000)

# Scaling the data and clustering
all.genes <- rownames(pbmc_sub)
pbmc_sub <- ScaleData(pbmc_sub) #features = all.genes)
pbmc_sub <- RunPCA(pbmc_sub, features = VariableFeatures(object = pbmc_sub))
pbmc_sub <- FindNeighbors(pbmc_sub, dims = 1:10)
pbmc_sub <- FindClusters(pbmc_sub, resolution = 0.5)

# Run UMAP/tSNE
pbmc_sub <- RunUMAP(pbmc_sub, dims = 1:10)
DimPlot(pbmc_sub, reduction = "umap")

# Save the file at this point
# saveRDS(pbmc_sub, file = "/Users/jamesbrenner/Desktop/rotations/Wu Lab/singe cell practice/filtered_gene_bc_matrices/hg19/pbmc_tutorial.rds")


## subset the T cells and re-cluster
pbmc_sub$random <- Idents(pbmc_sub)
Idents(pbmc_sub) <- "CELLTYPE"
pbmc_tcell <- subset(pbmc_sub, idents = c("CD4", "CD8", "Unconventional T"))
pbmc_tcell <- FindVariableFeatures(pbmc_tcell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell)
pbmc_tcell <- ScaleData(pbmc_tcell) 
pbmc_tcell <- RunPCA(pbmc_tcell, features = VariableFeatures(object = pbmc_tcell))
pbmc_tcell <- FindNeighbors(pbmc_tcell, dims = 1:10)
pbmc_tcell <- FindClusters(pbmc_tcell, resolution = 0.5)

# Run UMAP/tSNE
pbmc_tcell <- RunUMAP(pbmc_tcell, dims = 1:10)
DimPlot(pbmc_tcell, reduction = "umap", label = TRUE)

pbmc_tcell@meta.data

pbmc_tcell$CellType <- Idents(pbmc_tcell)
# Next, switch the identity class of all cells to reflect replicate ID
#Idents(pbmc_tcell) <- "CELLTYPE2"
Idents(pbmc_tcell) <- "seurat_clusters"
#Idents(pbmc_tcell) <- "Allpre_post"
#Idents(pbmc_tcell) <- "response"

DimPlot(pbmc_tcell, reduction = "umap")
ggsave(file = "plots/tcell_sd_umap_clusters.eps")

DotPlot(object = pbmc_tcell, features = c("CD4", "CD8A", "CCR7", "TIGIT", "GZMA", "GZMB", "ZNF683"), cols="RdBu") +
  theme(axis.text.x = element_text(angle=45)) +
  theme(axis.text.x = element_text(hjust=1))
ggsave(file = "plots/tcell_sd_dot_features.eps")


pbmc_tcell_c6 <- subset(pbmc_tcell, idents = "6")
pbmc_tcell_c6 <- FindNeighbors(pbmc_tcell_c6, dims = 1:10)
pbmc_tcell_c6 <- FindClusters(pbmc_tcell_c6, resolution = 0.5)
Idents(pbmc_tcell_c6) <- "CELLTYPE2"
DimPlot(pbmc_tcell_c6, reduction = "umap")
FeaturePlot(pbmc_tcell_c6, features = c("ZNF683"))


### Finding deferentially expressed genes across clusters 
pbmc.markers_tcell <- FindAllMarkers(pbmc_tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tcell_diff_exp <- pbmc.markers_tcell %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# Visualize expression markers across clusters 
VlnPlot(pbmc_tcell, features = c("HNRNPLL", "CD8A", "CCR7", "KLRG1", "CD27", "CXCR3", "CX3CR1", "CD4", "ZNF683")) #CD244
# you can plot raw counts as well
VlnPlot(pbmc_tcell, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# Plot expression over UMAPs
FeaturePlot(pbmc_tcell, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                                   "CD8A"))

pbmc.markers_tcell %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(pbmc_tcell, features = top10$gene) + NoLegend()

cluster2.markers <- FindMarkers(pbmc_tcell, ident.1 = 2, ident.2 = c(0, 1, 3, 4, 9), min.pct = 0.25)
head(cluster2.markers, n = 20)
cluster6.markers <- FindMarkers(pbmc_tcell, ident.1 = 6, ident.2 = c(0, 1, 3, 4, 9), min.pct = 0.25)
head(cluster6.markers, n = 20)
cluster7.markers <- FindMarkers(pbmc_tcell, ident.1 = 7, ident.2 = c(0, 1, 3, 4, 9), min.pct = 0.25)
head(cluster7.markers, n = 20)
cluster8.markers <- FindMarkers(pbmc_tcell, ident.1 = 8, ident.2 = c(0, 1, 3, 4, 9), min.pct = 0.25)
head(cluster8.markers, n = 20)

cluster2.markers_temra <- FindMarkers(pbmc_tcell, ident.1 = 2, ident.2 = c(6, 7, 8), min.pct = 0.25)
head(cluster2.markers_temra, n = 20)

### Determine number of cells recorded in each cluster 
table(Idents(pbmc_tcell), pbmc_tcell$response)
table(Idents(pbmc_tcell), pbmc_tcell$ident_nl)
numbers <- table(Idents(pbmc_tcell), pbmc_tcell$pre_post)

### Calculate frequency of cells in particular cluster
numbers2 <- data.frame(numbers)
numbers2 <- pivot_wider(numbers2, 
                        names_from = "Var1",
                        values_from = "Freq")
numbers2 <- numbers2 %>% mutate(sum = rowSums(across(where(is.numeric))))
numbers2$freq2 <- numbers2$"2" / numbers2$"sum"
numbers2$freq6 <- numbers2$"6" / numbers2$"sum"
numbers2$freq7 <- numbers2$"7" / numbers2$"sum"
write_csv(numbers2, file="plots/sd_numbers_by_cluster.csv")

to_graph <- select(numbers2, "Var2", "freq2")
to_graph6 <- select(numbers2, "Var2", "freq6")
to_graph7 <- select(numbers2, "Var2", "freq7")

r_pre <- c("PT1-Pre", "PT2-Pre", "PT3-Pre")
r_post <- c("PT1-Post", "PT2-Post", "PT3-Post")
nr_pre <- c("PT4-Pre", "PT5-Pre", "PT6-Pre", "PT7-Pre", "PT8-Pre")
nr_post <- c("PT4-Post", "PT5-Post", "PT6-Post", "PT7-Post", "PT8-Post")
nr_pre_nosd <- c("PT4-Pre", "PT5-Pre", "PT6-Pre")
nr_post_nosd <- c("PT4-Post", "PT5-Post", "PT6-Post")

to_graph <- to_graph[-1,]
to_graph <- separate(to_graph, "Var2", c("pat","pre_post"), "-", extra= "merge")
to_graph_r <- subset(to_graph, pat %in% c("PT1", "PT2", "PT3"))
to_graph_nr <- subset(to_graph, pat %in% c("PT4", "PT5", "PT6")) #, "PT7", "PT8"))
to_graph_sd <- subset(to_graph, pat %in% c("PT7", "PT8"))

to_graph6 <- to_graph6[-1,]
to_graph6 <- separate(to_graph6, "Var2", c("pat","pre_post"), "-", extra= "merge")
to_graph6_r <- subset(to_graph6, pat %in% c("PT1", "PT2", "PT3"))
to_graph6_nr <- subset(to_graph6, pat %in% c("PT4", "PT5", "PT6")) #, "PT7", "PT8"))
to_graph6_sd <- subset(to_graph6, pat %in% c("PT7", "PT8"))

to_graph7 <- to_graph7[-1,]
to_graph7 <- separate(to_graph7, "Var2", c("pat","pre_post"), "-", extra= "merge")
to_graph7_r <- subset(to_graph7, pat %in% c("PT1", "PT2", "PT3"))
to_graph7_nr <- subset(to_graph7, pat %in% c("PT4", "PT5", "PT6")) #, "PT7", "PT8"))
to_graph7_sd <- subset(to_graph7, pat %in% c("PT7", "PT8"))

to_graph_r <- to_graph_r %>%
  arrange(desc(pre_post), pat) %>%
  arrange((pat), pat)
to_graph_nr <- to_graph_nr %>%
  arrange(desc(pre_post), pat) %>%
  arrange((pat), pat)
to_graph_sd <- to_graph_sd %>%
  arrange(desc(pre_post), pat) %>%
  arrange((pat), pat)

### Graph frequencies 
to_graph_r %>%
  ggplot(aes(x=reorder(pre_post, freq2), y=freq2, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R", y = "Frequency of Cluster 2 Cells") +

  to_graph_nr %>%
  ggplot(aes(x=reorder(pre_post, -freq2), y=freq2, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "NR") + 
  
  to_graph_sd %>%
  ggplot(aes(x=reorder(pre_post, freq2), y=freq2, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "SD") 

to_graph6_r %>%
  ggplot(aes(x=reorder(pre_post, freq6), y=freq6, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R", y = "Frequency of Cluster 6 Cells") +
  
  to_graph6_nr %>%
  ggplot(aes(x=reorder(pre_post, freq6), y=freq6, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "NR") + 
  
  to_graph6_sd %>%
  ggplot(aes(x=reorder(pre_post, -freq6), y=freq6, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "SD") 

to_graph7_r %>%
  ggplot(aes(x=reorder(pre_post, freq7), y=freq7, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R", y = "Frequency of Cluster 7 Cells") +
  
  to_graph7_nr %>%
  ggplot(aes(x=reorder(pre_post, -freq7), y=freq7, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "NR") + 
  
  to_graph7_sd %>%
  ggplot(aes(x=reorder(pre_post, freq7), y=freq7, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "SD") 

prop.table(table(Idents(pbmc_tcell), pbmc_tcell$response), margin = 2)

pbmc_r <- subset(pbmc_tcell, subset = pre_post == r_pre)
plot8 <- VlnPlot(pbmc_r, features = c("HNRNPLL", "CD8A", "CCR7", "KLRG1", "CD27", "CXCR3", "CX3CR1", "CD4", "ZNF683")) #CD244
pbmc_nr <- subset(pbmc_tcell, subset = pre_post == r_post)
plot9 <- VlnPlot(pbmc_nr, features = c("HNRNPLL", "CD8A", "CCR7", "KLRG1", "CD27", "CXCR3", "CX3CR1", "CD4", "ZNF683")) #CD244
plot9

# Plot feature expression graphs over UMAP

FeaturePlot(pbmc_tcell, features = c("ZNF683"))

# Caluclate average expression per cluster for features 

combined_averages <- AverageExpression(pbmc_tcell, return.seurat = TRUE) 
r_pre_avg <- AverageExpression(subset(pbmc_tcell, subset = pre_post == r_pre), return.seurat = TRUE)
r_post_avg <- AverageExpression(subset(pbmc_tcell, subset = pre_post == r_post), return.seurat = TRUE)
nr_pre_avg <- AverageExpression(subset(pbmc_tcell, subset = pre_post == nr_pre), return.seurat = TRUE)
nr_post_avg <- AverageExpression(subset(pbmc_tcell, subset = pre_post == nr_post), return.seurat = TRUE)
r_avg <- AverageExpression(subset(pbmc_tcell, subset = pre_post == c(r_pre, r_post)), return.seurat = TRUE)
nr_avg <- AverageExpression(subset(pbmc_tcell, subset = pre_post == c(nr_pre, nr_post)), return.seurat = TRUE)

# Create heat map 

heat1 <- DoHeatmap(r_avg, features = c("HNRNPLL", "CD8A", "CCR7", "KLRG1", "CD27", "CXCR3", "CX3CR1", "CD4", "ZNF683"), 
            label = FALSE ,draw.lines = FALSE)  +  
            scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu")))
heat2 <- DoHeatmap(nr_avg, features = c("HNRNPLL", "CD8A", "CCR7", "KLRG1", "CD27", "CXCR3", "CX3CR1", "CD4", "ZNF683"), 
            label = TRUE ,draw.lines = FALSE)  +  
            scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu")))
#heat1 + heat2 

gene.list <- c("CD3E", "CD4", "CD8A", 
               "SELL", "IL7R", "CCR7", "CD27", "CD28", "HNRNPLL", "PTPRC",
               "CCL5", "FAS", "CD69", "ITGAE", "CXCR4", 
               "PDCD1", "KLRG1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", 
               "ENTPD1", "IL2RA", "ICOS", "TNFRSF4", "CD40LG", 
               "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "TNF", 
               "LEF1", "TCF7", "EOMES", "TBX21", "TOX", "ID2", "ID3", "ZNF683", "MKI67", 
               "TRGV9", "KLRB1", "KLRC1", "TRAV1-2")

num_idents <- length(unique(Idents(object = pbmc_tcell)))
na_colors <- rep("NA", num_idents)

heat3 <- DoHeatmap(nr_pre_avg, features = gene.list, label = TRUE ,draw.lines = FALSE, group.colors = na_colors) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "pre") 

heat4 <- DoHeatmap(nr_post_avg, features = gene.list, label = TRUE ,draw.lines = FALSE, group.colors = na_colors) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "post") 

heat3 + heat4 

heat14 <- DoHeatmap(combined_averages, features = gene.list, label = TRUE ,draw.lines = FALSE, group.colors = na_colors) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "post") 
heat14


### subset without SD patients 
pbmc_tcell_test <- subset(pbmc_tcell, subset = response == c("CRPR", "NR"))
pbmc_tcell_test <- FindVariableFeatures(pbmc_tcell_test, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell_test)
pbmc_tcell_test <- ScaleData(pbmc_tcell_test) 
pbmc_tcell_test <- RunPCA(pbmc_tcell_test, features = VariableFeatures(object = pbmc_tcell_test))
pbmc_tcell_test <- FindNeighbors(pbmc_tcell_test, dims = 1:10)
pbmc_tcell_test <- FindClusters(pbmc_tcell_test, resolution = 0.5)

# Run UMAP/tSNE
pbmc_tcell_test <- RunUMAP(pbmc_tcell_test, dims = 1:10)
DimPlot(pbmc_tcell_test, reduction = "umap", label = TRUE)
ggsave(file = "plots/tcell_nosd_umap_clusters.eps", width = 8, height = 6, dpi = 150)
DimPlot(pbmc_tcell_test, reduction = "umap", label = FALSE)
ggsave(file = "plots/tcell_nosd_umap_clusters_nolab.eps", width = 8, height = 6, dpi = 150)

### Save file 
saveRDS(pbmc_tcell_test, file = "objects/nosd_tcells.rds")

Idents(pbmc_tcell_test) <- "CELLTYPE2"
DimPlot(pbmc_tcell_test, reduction = "umap", label = FALSE)
ggsave(file = "plots/tcell_nosd_umap_celltype_nolab.eps", width = 8, height = 6, dpi = 150)

Idents(pbmc_tcell_test) <- "response"
DimPlot(pbmc_tcell_test, reduction = "umap", label = FALSE)
ggsave(file = "plots/tcell_nosd_umap_response_nolab.eps", width = 8, height = 6, dpi = 150)

Idents(pbmc_tcell_test) <- "Allpre_post"
DimPlot(pbmc_tcell_test, reduction = "umap", label = FALSE)
ggsave(file = "plots/tcell_nosd_umap_pre_post_nolab.eps", width = 8, height = 6, dpi = 150)

Idents(pbmc_tcell_test) <- "seurat_clusters"

### Adding TCR data
tcrs <- read.csv("TCRs_AML_Abbas.csv", header = TRUE)
meta3 <- pbmc_tcell_test@meta.data
pbmc_tcell_test@meta.data <- merge(pbmc_tcell_test@meta.data, tcrs, by.x = "V1", by.y = "cellID",all.x = TRUE, sort=FALSE)
rownames(pbmc_tcell_test@meta.data) <- pbmc_tcell_test@meta.data$V1
pbmc_tcell_test@meta.data <- pbmc_tcell_test@meta.data[meta3$V1,]
rownames(pbmc_tcell_test@meta.data) <- rownames(meta3)

head(pbmc_tcell_test@meta.data)
Idents(pbmc_tcell_test) <- "seurat_clusters"
DimPlot(pbmc_tcell_test, reduction = "umap", label = TRUE)
FeaturePlot(pbmc_tcell_test, features = "clone.size")
FeaturePlot(pbmc_tcell_test, features = "BCL11B")
FeaturePlot(pbmc_tcell_test, features = "KLRC2") #for NKG2C
FeaturePlot(pbmc_tcell_test, features = "KLRC1")

pbmc_tcell_test@meta.data$log2 <- log2(pbmc_tcell_test@meta.data$clone.size)
FeaturePlot(pbmc_tcell_test, features = "log2")
ggsave(file = "plots/tcell_nosd_umap_clonal_expansion.eps", width = 8, height = 6, dpi = 150)


##### Trying to make a pie graph, unsuccessful, please disregard 
pbmc_tcell_test@meta.data
pies <- subset(pbmc_tcell_test, seurat_clusters == 2 & response == "CRPR" & Allpre_post == "Pre") 
pies2 <- subset(pbmc_tcell_test, seurat_clusters == 2 & response == "CRPR" & Allpre_post == "Post") 
pies2@meta.data

pie_numbers <- as.data.frame(pies@meta.data$patient.clonotypeID)
colnames(pie_numbers)[1] ="values"
pie_numbers$ID <- "pre"
pie2_numbers <- as.data.frame(pies2@meta.data$patient.clonotypeID)
colnames(pie2_numbers)[1] ="values"
pie2_numbers$ID <- "post"
pie_all <- bind_rows(pie_numbers, pie2_numbers)
pie_all <- drop_na(pie_all)
pie_all$number <- 1


pie_all_num <- pie_all %>%
  group_by(values, ID) %>%
  summarise(sum = sum(number))
write_csv(pie_all_num, file="plots/clonotype_expansion_numbers.csv")


pie_all_exp <- pie_all_num[sum > 1 & ID == "post"] 

pie_all_num %>%
  ggplot(aes(x=reorder(ID, sum), y=sum, fill=ID)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=values)) +
  scale_y_continuous(limits = c(0, 40)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "pie") 

pie(table(pie_numbers$values)) + 
  theme_void()
  
pie(table(pie2_numbers$values))

ggplot() +
  geom_bar(data = pie_all_exp, aes(x=ID, y=sum, fill=values), width = 0.2, stat="identity") +
  theme(legend.position = "none") 

  geom_point()+ 
#  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R_test", y = "Frequency of Cluster 2 Cells")

  
### Extract data for frequency calculations 
numbers_test <- table(Idents(pbmc_tcell_test), pbmc_tcell_test$pre_post)

FeaturePlot(pbmc_tcell_test, features = c("CD8A", "CD4", "ZNF683", "KLRC2"))
ggsave(file = "plots/tcell_nosd_umap_feature_expression.eps", width = 8, height = 6, dpi = 150)


numbers_test_bytime <- table(Idents(pbmc_tcell_test), pbmc_tcell_test$orig.ident)
yes <- data.frame(numbers_test_bytime)
yes <- pivot_wider(yes, 
                             names_from = "Var1",
                             values_from = "Freq")


num_idents_test <- length(unique(Idents(object = pbmc_tcell_test)))
na_colors_test <- rep("NA", num_idents_test)

combined_averages_test <- AverageExpression(pbmc_tcell_test, return.seurat = TRUE) 
r_pre_avg_test <- AverageExpression(subset(pbmc_tcell_test, subset = pre_post == r_pre), return.seurat = TRUE)
r_post_avg_test <- AverageExpression(subset(pbmc_tcell_test, subset = pre_post == r_post), return.seurat = TRUE)
nr_pre_nosd_avg_test <- AverageExpression(subset(pbmc_tcell_test, subset = pre_post == nr_pre), return.seurat = TRUE)
nr_post_nosd_avg_test <- AverageExpression(subset(pbmc_tcell_test, subset = pre_post == nr_post), return.seurat = TRUE)
r_avg_test <- AverageExpression(subset(pbmc_tcell_test, subset = pre_post == c(r_pre, r_post)), return.seurat = TRUE)
nr_avg_test <- AverageExpression(subset(pbmc_tcell_test, subset = pre_post == c(nr_pre, nr_post)), return.seurat = TRUE)


heat5 <- DoHeatmap(r_avg_test, features = gene.list, label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "R") 

heat6 <- DoHeatmap(nr_avg_test, features = gene.list, label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "NR") 

heat5 + heat6

### Calculate the frequencies of clusters 
numbers2_test <- data.frame(numbers_test)
numbers2_test <- pivot_wider(numbers2_test, 
                        names_from = "Var1",
                        values_from = "Freq")
numbers2_test <- numbers2_test %>% mutate(sum = rowSums(across(where(is.numeric))))
numbers2_test$freq2 <- numbers2_test$"2" / numbers2_test$"sum"
numbers2_test$freq6 <- numbers2_test$"6" / numbers2_test$"sum"
numbers2_test$freq7 <- numbers2_test$"7" / numbers2_test$"sum"
write_csv(numbers2_test, file="plots/numbers_by_cluster_nosd.csv")
numbers2_test <- read_csv(file = '/Users/jamesbrenner/My Drive/final_data/aml_pd1/plots/numbers_by_cluster_nosd.csv')

to_graph_test <- select(numbers2_test, "Var2", "freq2")
to_graph6_test <- select(numbers2_test, "Var2", "freq6")
to_graph7_test <- select(numbers2_test, "Var2", "freq7")

r_pre <- c("PT1-Pre", "PT2-Pre", "PT3-Pre")
r_post <- c("PT1-Post", "PT2-Post", "PT3-Post")
nr_pre <- c("PT4-Pre", "PT5-Pre", "PT6-Pre", "PT7-Pre", "PT8-Pre")
nr_post <- c("PT4-Post", "PT5-Post", "PT6-Post", "PT7-Post", "PT8-Post")
nr_pre_nosd <- c("PT4-Pre", "PT5-Pre", "PT6-Pre")
nr_post_nosd <- c("PT4-Post", "PT5-Post", "PT6-Post")

to_graph_test <- separate(to_graph_test, "Var2", c("pat","pre_post"), "-", extra= "merge")
to_graph_r_test <- subset(to_graph_test, pat %in% c("PT1", "PT2", "PT3"))
to_graph_nr_test <- subset(to_graph_test, pat %in% c("PT4", "PT5", "PT6")) 

to_graph6_test <- separate(to_graph6_test, "Var2", c("pat","pre_post"), "-", extra= "merge")
to_graph6_r_test <- subset(to_graph6_test, pat %in% c("PT1", "PT2", "PT3"))
to_graph6_nr_test <- subset(to_graph6_test, pat %in% c("PT4", "PT5", "PT6"))

to_graph7_test <- separate(to_graph7_test, "Var2", c("pat","pre_post"), "-", extra= "merge")
to_graph7_r_test <- subset(to_graph7_test, pat %in% c("PT1", "PT2", "PT3"))
to_graph7_nr_test <- subset(to_graph7_test, pat %in% c("PT4", "PT5", "PT6"))

to_graph_r_test <- to_graph_r_test %>%
  arrange(desc(pre_post), pat) %>%
  arrange((pat), pat)
to_graph_nr_test <- to_graph_nr_test %>%
  arrange(desc(pre_post), pat) %>%
  arrange((pat), pat)


### Graph frequencies 
to_graph_r_test %>%
  ggplot(aes(x=reorder(pre_post, freq2), y=freq2, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R", y = "Frequency of Cluster 2 Cells") +
  stat_compare_means(method = 't.test', 
                     paired = TRUE, 
                     comparisons = list( c('Pre', 'Post'))) +
  
  to_graph_nr_test %>%
  ggplot(aes(x=reorder(pre_post, freq2), y=freq2, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "NR") +
  stat_compare_means(method = 't.test', 
                     paired = TRUE, 
                     comparisons = list( c('Pre', 'Post')))
ggsave(file = "plots/tcell_nosd_cluster2_bar.eps", width = 5, height = 6, dpi = 150)

to_graph6_r_test %>%
  ggplot(aes(x=reorder(pre_post, -freq6), y=freq6, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R", y = "Frequency of Cluster 6 Cells") +
  
  to_graph6_nr_test %>%
  ggplot(aes(x=reorder(pre_post, -freq6), y=freq6, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "NR") 
ggsave(file = "plots/tcell_nosd_cluster6_bar.eps", width = 5, height = 6, dpi = 150)

to_graph7_r_test %>%
  ggplot(aes(x=reorder(pre_post, freq7), y=freq7, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  labs(title = "R", y = "Frequency of Cluster 7 Cells") +
  
  to_graph7_nr_test %>%
  ggplot(aes(x=reorder(pre_post, -freq7), y=freq7, fill=pre_post)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=pat)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "NR") 
ggsave(file = "plots/tcell_nosd_cluster7_bar.eps", width = 5, height = 6, dpi = 150)
  
### Dot and heat plots of features 
DotPlot(object = pbmc_tcell_test, features = c("CD4", "CD8A", "CCR7", "TIGIT", "GZMA", "GZMB", "ZNF683"), cols="RdBu") +
  theme(axis.text.x = element_text(angle=45)) +
  theme(axis.text.x = element_text(hjust=1))
dot1 <- DotPlot(object = subset(pbmc_tcell_test, subset = pre_post == c(r_pre, r_post)), features = c("CD4", "CD8A", "CCR7", "TIGIT", "GZMA", "GZMB", "GZMK", "ZNF683"), cols="RdBu") +
  theme(axis.text.x = element_text(angle=45)) +
  theme(axis.text.x = element_text(hjust=1)) +
  labs(title = "R") 
dot2 <- DotPlot(object = subset(pbmc_tcell_test, subset = pre_post == c(nr_pre, nr_post)), features = c("CD4", "CD8A", "CCR7", "TIGIT", "GZMA", "GZMB", "GZMK", "ZNF683"), cols="RdBu") +
  theme(axis.text.x = element_text(angle=45)) +
  theme(axis.text.x = element_text(hjust=1)) +
  labs(title = "NR") 
dot1 + dot2 
ggsave(file = "plots/tcell_nosd_dot_rnr.eps", width = 10, height = 6, dpi = 150)

heat7 <- DoHeatmap(combined_averages_test, features = gene.list, label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "Clustered") 
heat7 
ggsave(file = "plots/tcell_nosd_heat.eps", width = 6, height = 10, dpi = 150)


### subset with NK cells 
Idents(pbmc_sub) <- "CELLTYPE"
pbmc_tcell_nk <- subset(pbmc_sub, idents = c("CD4", "CD8", "Unconventional T", "preT/NK", "NK"))
pbmc_tcell_nk <- subset(pbmc_tcell_nk, subset = response == c("CRPR", "NR"))
pbmc_tcell_nk <- FindVariableFeatures(pbmc_tcell_nk, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_tcell_test)
pbmc_tcell_test <- ScaleData(pbmc_tcell_test) 
pbmc_tcell_test <- RunPCA(pbmc_tcell_test, features = VariableFeatures(object = pbmc_tcell_test))
pbmc_tcell_nk <- FindNeighbors(pbmc_tcell_nk, dims = 1:10)
pbmc_tcell_nk <- FindClusters(pbmc_tcell_nk, resolution = 1)

# Run UMAP/tSNE across various descriptors 
pbmc_tcell_nk <- RunUMAP(pbmc_tcell_nk, dims = 1:10)
Idents(pbmc_tcell_nk) <- "seurat_clusters"
DimPlot(pbmc_tcell_nk, reduction = "umap", label = TRUE)
ggsave(file = "plots/tcell_nk_umap_clusters.eps", width = 8, height = 6, dpi = 150)

### Save file 
saveRDS(pbmc_tcell_nk, file = "objects/nosd_tcells_nks.rds")

Idents(pbmc_tcell_nk) <- "CELLTYPE2"
DimPlot(pbmc_tcell_nk, reduction = "umap", label = FALSE)
ggsave(file = "plots/tcell_nk_umap_celltype_nolab.eps", width = 8, height = 6, dpi = 150)

FeaturePlot(pbmc_tcell_nk, features = c("ZNF683", "CD8A", "CD4", "IL2RB"))
ggsave(file = "plots/tcell_nk_umap_features.eps", width = 8, height = 6, dpi = 150)
Idents(pbmc_tcell_nk) <- "seurat_clusters"

#Idents(pbmc_tcell_nk) <- "response"
#Idents(pbmc_tcell_nk) <- "orig.ident"
#Idents(pbmc_tcell_nk) <- "Allpre_post"


### Dot and heat plots of features 
combined_averages_nk <- AverageExpression(pbmc_tcell_nk, return.seurat = TRUE) 
gene.list.nk <- c("ZNF683", "CD8A", "CD4", "CD3E", "BCL11B", 
                  "SYK", "KIR2DL3", "NCR1", 
                  "GZMB", "KIR2DL1", "KIR2DL4", "KIR3DL2",
                  "KLRC2", "KLRC3", "TYROBP", "CD244", 
                  "KLRC1", "KIR3DL1", "PLCG2", "FCGR3A",
                  "ITGAX", "S1PR5", "PRF1", "TBX21", "GNLY",
                  "IL2RB", "LAT2", "LYN", "NCAM1", "SLAMF7", 
                  "TCF7", "CD27", "CD5", "PDCD1", "CD28", 
                  "CCR7", "CD6", "IL7R", "THEMIS", "TESPA1", 
                  "B3GAT1", "XCL1", "XCL2", "CD160", "SELL")

gene.list.tcell.cd8 <- c("CD8A", "CD4", "CD69", "ITGAE", 
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

DoHeatmap(combined_averages_nk, features = gene.list.nk, label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "Clustered") 
ggsave(file = "plots/tcell_nk_heat.eps", width = 6, height = 10, dpi = 150)

DotPlot(object = pbmc_tcell_nk, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/tcell_nk_dot_cd8_panel.eps", width = 15, height = 8, dpi = 150,)

DotPlot(object = pbmc_tcell_test, features = gene.list.tcell.cd8, cols="RdBu") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.text.x = element_text(hjust=1)) +
  theme(axis.text.x = element_text(vjust=0.5)) +
  theme(axis.text.y = element_text(angle=90)) +
  theme(axis.text.y = element_text(hjust=0.5)) 
ggsave(file = "plots/tcell_nosd_dot_cd8_panel.eps", width = 15, height = 8, dpi = 150,)

cluster2.markers.nk <- FindMarkers(pbmc_tcell_nk, ident.1 = 2, min.pct = 0.25)#ident.2 = c(3, 4), 
head(cluster2.markers.nk, n = 20)

DoHeatmap(combined_averages_test, features = gene.list.nk, label = TRUE ,draw.lines = FALSE, group.colors = na_colors_test) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =8, name = "RdBu"))) +
  labs(title = "Clustered") 






