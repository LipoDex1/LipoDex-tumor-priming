# Publicly available code for "Desmoplastic Tumor Priming using Clinical-stage Dexamethasone Liposomes"
# Written by Gideon JL Schaefer, RWTH Aachen University

# ============================== Quality Control ============================== #

# Load required libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(ggpubr)
library(writexl)
library(openxlsx)
library(dittoSeq)
library(stringr)
library(magick)
library(svglite)

# Define object and input path
sample <- "sample_name"
input_path <- "data"  # Path to /filtered_feature_bc_matrix/

# Import dataset
data <- Read10X(data.dir = input_path)
data <- CreateSeuratObject(counts = data, project = "FP6", min.cells = 3, min.features = 200)

# Compute mitochondrial and ribosomal RNA percentages
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = '^mt-', assay = "RNA")
data[["percent.rib"]] <- PercentageFeatureSet(data, pattern = '^Rpl|^Rps', assay = "RNA")

# Plot distributions
VlnPlot(data, features = "percent.rib")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)

# Cell filtering criteria
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & percent.rib < 50)

# Add metadata
data@meta.data$treatment <- "treatment"  # Lipodexamethasone vs. Control

# ============================== Harmony Integration ============================== #

# Load required libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(clustree)
library(genesorteR, quietly = TRUE)
library(dplyr)
library(patchwork)
library(harmony)

# Load sample data
sample <- readRDS(file = "samplepath.RDS")

# Merge samples, normalize, and scale data
merged <- merge(sample1, y = c(sample2), add.cell.ids = c("sample1", "sample2"), project = "merged")
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)

# Run Harmony integration
merged <- RunHarmony(merged, "treatment", plot_convergence = TRUE)

# ============================== Cluster Annotation ============================== #

# Load required libraries
library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(clustree)
library(dplyr)
library(patchwork)
library(harmony)

# Load merged Harmony object
data <- readRDS(file = "merged")

# Define cluster identities
Idents(data) <- factor(data$treatment, levels = c("Control", "LipoDex"))

# Find variable features
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 4000)

# Perform clustering with predefined parameters
data <- NormalizeData(data)
data <- ScaleData(data)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:25)
data <- FindClusters(data, resolution = 0.1)
data <- RunUMAP(data, reduction = "harmony", dims = 1:25)

# UMAP visualization
DimPlot(data, reduction = "umap", label = FALSE, label.size = 10, pt.size = 0.1)

# Identify cluster markers
all.markers <- FindAllMarkers(data, only.pos = TRUE)
all.markers <- all.markers %>% group_by(seurat_clusters) %>% dplyr::filter(avg_log2FC > 1)

# Rename cluster identities
data <- RenameIdents(data, "0" = "MPC", "1" = "TAM", "2" = "prolif_MPC", "3" = "fibroblasts", "4" = "epithelial", "5" = "endothelial")
data$clusternames <- Idents(data)

# Generate dot plot with top genes
clusters <- 0:5
top5_all <- do.call(rbind, lapply(clusters, function(c) {
  head(all.markers[all.markers$cluster == c,], 5)
}))
Idents(data) <- factor(data$clusternames, levels = c("MPC", "TAM", "prolif_MPC", "fibroblasts", "epithelial", "endothelial"))
DotPlot(data, features = top5_all$gene) + RotatedAxis() + coord_flip()

# Final UMAP
DimPlot(data, label = TRUE, label.size = 5, pt.size = 0.1) + NoLegend()

# ============================== Cluster Annotation ============================== #
# Load required libraries
library(RColorBrewer)
library(ggpubr)
library(clustree)
library(genesorteR, quietly = TRUE)
library(dplyr)
library(patchwork)
library(harmony)
library(dittoSeq)
library(slingshot)
library(tidymodels)
library(writexl)
library(DESeq2)
library(enrichR)
library(openxlsx)
library(ggrepel)
library(EnhancedVolcano)
library(AnnotationDbi)
library(gprofiler2)
BiocManager::install("org.Mm.eg.db")

# Define object
data <- readRDS(file = "data.RDS") ## Pass cluster annotated harmony integrated Seurat object
data@meta.data$clusternames <- Idents(data)

# Composition visualization
dittoBarPlot(data, var = Idents(data), group.by = "treatment", scale = "percent", data.out = FALSE, ylab = "Ratio of Cellcluster", color.panel = c("darkblue", "yellow", "red", "skyblue3", "darkorange", "darkgreen"))
Idents(data) <- data$clusternames
data <- subset(data, idents = c("prolif_MPC", "MPC", "fibroblasts", "TAM"))
dittoBarPlot(data, var = Idents(data), group.by = "treatment", scale = "percent", data.out = FALSE, ylab = "Cluster Ratio")

# Find DEGs; definition: upregulated (avg_log2FC > 0) / downregulated (avg_log2FC < 0) in treatment group (ident.1 = LipoDex) compared to control group (ident.2 = Control)
data$celltype_treatment <- paste0(data$clusternames, "_", data$treatment)
Idents(data) <- data$celltype_treatment

#______TAM______#
TAM_de <- FindMarkers(data, ident.1 = "TAM_LipoDex", ident.2 = "TAM_Control", verbose = FALSE, only.pos = FALSE, test.use = "wilcox")
TAM_de$gene <- rownames(TAM_de)
TAM_de <- TAM_de[TAM_de$p_val_adj <= 0.05,] 
GOI_TAM_up <- c("Ctu2", "Stk4", "Pld1", "Nrf1", "Cyld", "Cse1l", "Klf10", "Otulin", "Ptpn9", "Tnip1")
GOI_TAM_down <- c("Parp12", "Sos1", "Tradd", "Nfkb2", "Egln2", "Cdk9", "Akt1", "Micu2")
EnhancedVolcano(TAM_de, 
                x = "avg_log2FC", 
                y = "p_val_adj", 
                lab = TAM_de$gene, 
                FCcutoff = 0.15, 
                pCutoff = 0.05, 
                xlab = "Fold Change (log2)", 
                pointSize = 0.8,  
                labSize = 6, 
                xlim = c(-0.5, 2.5), 
                ylim = c(0, 11), 
                selectLab = "none", 
                legendPosition = "none") +
  geom_text_repel(
    aes(label = ifelse(TAM_de$gene %in% c(GOI_TAM_up, GOI_TAM_down), TAM_de$gene, "")),
    box.padding = 0.5, 
    point.padding = 0.5, 
    segment.color = "grey50", 
    max.overlaps = Inf, 
    size = 6
  )

#______MPC______#
data_mpc <- data
data_mpc[["RNA"]]@layers$counts <- as.matrix(data[["RNA"]]@layers$counts)+1
MPC_de <- FindMarkers(data_mpc, ident.1 = "MPC_LipoDex", ident.2 =  "MPC_Control", only.pos = F, test.use = "wilcox")
MPC_de$gene <- rownames(MPC_de)
MPC_de <- MPC_de[MPC_de$p_val_adj<=0.05,] 

MPC_upreg_GOI <- c("Mmp3","Mrc2","Ptger2","Mmp2","F7","Lag3","Klf10","Klf11","Bap1")
MPC_downreg_GOI <- c("Nlrc3","Csf1","Ifi44","Angpt2","Ptges2","Tk1","Ptges","Ccl24","Csrp2","Lrrk2","Ddah2", "Ccr6", "Jak3", "Kdr", "Bgn", "Ecm1", "Fibp", "Fgf9", "Col1a2")


EnhancedVolcano(MPC_de, x = "avg_log2FC", y =  "p_val_adj", lab = MPC_de$gene, FCcutoff = 0.15, 
                pCutoff = 0.05, xlab = "Fold Change (log2)", pointSize = 0.8, labSize = 6, 
                xlim = c(-0.5, 1), ylim = c(0,120), selectLab = "none", 
                legendPosition = "none") +
  geom_text_repel(aes(label = ifelse(MPC_de$gene %in% c(MPC_upreg_GOI, MPC_downreg_GOI), MPC_de$gene, "")),
                  box.padding = 0.5, point.padding = 0.5, segment.color = "grey50", max.overlaps = Inf,size = 4)

#______prolif_MPC______#
prolif_MPC_de <- FindMarkers(data, ident.1 = "prolif_MPC_LipoDex", ident.2 =  "prolif_MPC_Control", verbose = FALSE, only.pos = F, test.use = "wilcox")
prolif_MPC_de$gene <- rownames(prolif_MPC_de)
prolif_MPC_de <- prolif_MPC_de[prolif_MPC_de$p_val_adj<=0.05,] 

#______fibroblasts______#
fibro_de <- FindMarkers(data, ident.1 = "fibroblasts_LipoDex", ident.2 = "fibroblasts_Control", verbose = FALSE, only.pos = F, test.use = "wilcox")
fibro_de$gene <- rownames(fibro_de)
fibro_de <- fibro_de[fibro_de$p_val_adj<=0.05,] 
Fibro_genes <- c("Fn1", "Col6a3", "Serpinh1", "Col1a1", "Lox", "Col5a2", "Col5a1", "Col6a1", "Col1a2", "Fbn1", "Col3a1", "Col12a1", "Fbln1", "Tnc", "Cd44", "Vcan", "Cxcl1", "Csf1")
EnhancedVolcano(fibro_de, x = "avg_log2FC", y =  "p_val_adj", lab = fibro_de$gene, FCcutoff = 0.05, pCutoff = 0.05, xlab = "Fold Change (log2)", pointSize = 0.5, labSize = 4, xlim = c(-2.2, 2), ylim = c(0,10)) +
  geom_text_repel(aes(label = ifelse(fibro_de$gene %in% Fibro_genes, fibro_de$gene, "")),
                  box.padding = 0.5, point.padding = 0.5, segment.color = "grey50", max.overlaps = Inf,size = 4)


# Gene Set enrichtment analysis with g Profiler
mpc_up <- rownames(MPC_de[MPC_de$avg_log2FC>0.1,])
mpc_down <- rownames(MPC_de[MPC_de$avg_log2FC<(-0.1),])
prolif_mpc_up <- rownames(prolif_MPC_de[prolif_MPC_de$avg_log2FC>0.1,])   
prolif_mpc_down <- rownames(prolif_MPC_de[prolif_MPC_de$avg_log2FC<(-0.1),])   
TAM_up <- rownames(TAM_de[TAM_de$avg_log2FC>0.1,])   
TAM_down <- rownames(TAM_de[TAM_de$avg_log2FC<(-0.1),])   
fibro_up <- rownames(fibro_de[fibro_de$avg_log2FC>0.1,])  
fibro_down <- rownames(fibro_de[fibro_de$avg_log2FC<(-0.1),])  

###MPC UP
GO_list <- gost(query = mpc_up, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_up_MPC <- result_GO
###MPC Down
GO_list <- gost(query = mpc_down, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_down_MPC <- result_GO
###prolif_MPC up
GO_list <- gost(query = prolif_mpc_up, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_up_prolif_MPC <- result_GO
result_GO <- shorten_strings(result_GO, intersection, max_length)
result_GO <- shorten_strings(result_GO, evidence_codes, max_length)
###prolif_MPC down
GO_list <- gost(query = prolif_mpc_down, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_down_prolif_MPC <- result_GO
result_GO <- shorten_strings(result_GO, intersection, max_length)
result_GO <- shorten_strings(result_GO, evidence_codes, max_length)
GO_list <- gost(query = TAM_up, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_up_TAM <- result_GO
result_GO <- shorten_strings(result_GO, intersection, max_length)
result_GO <- shorten_strings(result_GO, evidence_codes, max_length)
###TAM down
GO_list <- gost(query = TAM_down, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_down_TAM <- result_GO
result_GO <- shorten_strings(result_GO, intersection, max_length)
result_GO <- shorten_strings(result_GO, evidence_codes, max_length)
###fibro up
GO_list <- gost(query = fibro_up, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_up_fibro <- result_GO
result_GO <- shorten_strings(result_GO, intersection, max_length)
result_GO <- shorten_strings(result_GO, evidence_codes, max_length)
###fibro down
GO_list <- gost(query = fibro_down, organism = "mmusculus", ordered_query = F, significant = T, measure_underrepresentation = F, user_threshold = 0.05, evcodes = T)  
result_GO <- as.data.frame(GO_list$result)
GO_down_fibro <- result_GO
result_GO <- shorten_strings(result_GO, intersection, max_length)
result_GO <- shorten_strings(result_GO, evidence_codes, max_length)

#extract GO::BP 
GO_up_MPC_subset <- GO_up_MPC %>% 
  filter(source %in% c("GO:BP"))
GO_down_MPC_subset <- GO_down_MPC %>% 
  filter(source %in% c("GO:BP"))
GO_up_prolif_MPC_subset <- GO_up_prolif_MPC %>% 
  filter(source %in% c("GO:BP"))
GO_down_prolif_MPC_subset <- GO_down_prolif_MPC %>% 
  filter(source %in% c("GO:BP"))
GO_up_TAM_subset <- GO_up_TAM %>% 
  filter(source %in% c("GO:BP"))
GO_down_TAM_subset <- GO_down_TAM %>% 
  filter(source %in% c("GO:BP"))
GO_up_fibro_subset <- GO_up_fibro %>% 
  filter(source %in% c("GO:BP"))
GO_down_fibro_subset <- GO_down_fibro %>% 
  filter(source %in% c("GO:BP"))

