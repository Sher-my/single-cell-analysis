###Sample:Funiculus_umbilicalis
library(Seurat)
library(dplyr)
library(cowplot)
library(scran)
library(monocle) 

file_path <- "/home/lf-3/yt/rawdatas/sc_10/"
file_total <- list.files(file_path, pattern = "*.csv", full.names = TRUE)
samplename <- str_split_fixed(basename(file_total), "\\.c", n = 2)[,1] 
for (i in length(samplename)){}

BD_HC = read.delim2("/home/D:/My_data/single_cell/sc/K_FKDL202552766-1a_1_SampleTag02_hs_umbilical_cord_RSEC_MolsPerCell.csv", row.names = 1, header = TRUE, skip = 8, stringsAsFactors = FALSE, check.names = FALSE, sep = ",") 
BD_HC = t(BD_HC)
pbmc <- Seurat::CreateSeuratObject(counts = BD_HC, project = "umbilical_cord", min.cells = 3, min.features = 200)

##QC
pbmc[["percent.MT"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
HB.ref <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- intersect(HB.ref, rownames(pbmc@assays$RNA))
pbmc[["percent.HB"]] <- Seurat::PercentageFeatureSet(pbmc, features = HB.genes)
pbmc[["percent.Ribosome"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
plot1_1 <- Seurat::VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.HB", "percent.Ribosome"), ncol = 5)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/pre_qc.pdf', plot1_1, width = 12, height = 5, limitsize = FALSE)

##Normalize
qc_nFeature_RNA <- mean(pbmc$nFeature_RNA) + 3*sd(pbmc$nFeature_RNA)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < qc_nFeature_RNA & percent.MT < 35 & percent.HB < 0.2)
pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
plot2 <- Seurat::FeatureScatter(taipan_pbmc, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot3_1 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/c&p.png', plot2, width = 350, height = 300, limitsize = FALSE, units = ('mm'))
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/c&f.pdf', plot3_1, width = 5, height = 2, limitsize = FALSE)
plot2_1 <- Seurat::VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.HB", "percent.Ribosome"), ncol = 5)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/qc.pdf', plot2_1, width = 12, height = 5, limitsize = FALSE)

##Find var gene
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(Seurat::VariableFeatures(pbmc), 10)
plot4_1 <- Seurat::VariableFeaturePlot(pbmc)
plot5_1 <- Seurat::LabelPoints(plot = plot4_1, points = top10, repel = TRUE)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/all_g.pdf', plot4_1, width = 4, height = 3, limitsize = FALSE)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/top10_g.pdf', plot5_1, width = 6, height = 3, limitsize = FALSE)

##Cell cycle
all.genes <- rownames(pbmc)
pbmc <- Seurat::ScaleData(pbmc, features = all.genes, var.to.regress = "percent.MT")
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes
pbmc1 <- Seurat::CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc1 <- Seurat::RunPCA(pbmc1, features = c(s.genes, g2m.genes))
plot6_1 <- Seurat::DimPlot(pbmc1)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/cc.pdf', plot6_1, width = 5, height = 3, limitsize = FALSE)

pbmc1 <- Seurat::ScaleData(pbmc1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc1))
pbmc1 <- Seurat::RunPCA(pbmc1, features = c(s.genes, g2m.genes))
plot7 <- Seurat::DimPlot(pbmc1)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/decc.pdf', plot7, width = 5, height = 3, limitsize = FALSE)

##Pbmc
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
plot8_1 <- Seurat::ElbowPlot(pbmc, ndims = 50)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/dew.pdf', plot8_1, width = 4, height = 2, limitsize = FALSE)

##Cluster
pbmc <- Seurat::FindNeighbors(pbmc)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:15)
plot9_1 <- Seurat::DimPlot(pbmc, reduction = "umap")
taipan_pbmc <- RunTSNE(taipan_pbmc, dims = 1:30)
plot10 <- DimPlot(taipan_pbmc, reduction = "tsne")
saveRDS(taipan_pbmc, file = "/home/lf-3/yt/rawdatas/sc_10/pbmc.rds")
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/con-divide.pdf', plot9_1, width = 5, height = 3, limitsize = FALSE)

##Screen the differential expression genes
new.cluster.ids <- c("F","VCT", "CTB_cells", "EVT", "T-cells", "Macrophages", "DC", "MO")
names(new.cluster.ids) <- levels(taipan_pbmc)
taipan_pbmc <- Renameldents(pbmc, new.cluster.ids)
plot14 <- DimPlot(taipan_pbmc, reduction = "umap", label = TRUE, Pt.size = 0.5) + NoLegend()
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/sub&name.png', plot15, width = 350, height = 300, limitsize = FALSE, units = ('mm'))

pbmc_sub <- subset(x = pbmc, subset = CD47_isoform1 > 1)
plot15 <- Seurat::DimPlot(pbmc_sub, reduction = "umap", label = T)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/CD47_isoform1.pdf', plot15, width = 4, height = 3, limitsize = FALSE)

pbmc.markers <- Seurat::FindAllMarkers(object = pbmc, logfc.threshold = 0.25, min.pct = 0.25)
plot11 <- Seurat::VlnPlot(object = pbmc, features = c('CD47_isoform1'))
plot12 <- Seurat::FeaturePlot(object = pbmc, features = c('CD47_isoform1'))
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/CD47_isoform1_V.pdf', plot11, width = 4, height = 3, limitsize = FALSE)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/CD47_isoform1_F.pdf', plot12, width = 4, height = 3, limitsize = FALSE)

top10 <- taipan_pbmc.markers%>%group_by(cluster)%>%top_n(n = 10, wt = avg_logFC)
top10 <- as.data.frame(top10)
plot13 <- DoHeatmap(object = taipan_pbmc, feature = top10$gene)
ggplot2::ggsave('/home/lf-3/yt/rawdatas/sc_10/var_g.pdf', plot13, width = 1, height = 20, limitsize = FALSE)


