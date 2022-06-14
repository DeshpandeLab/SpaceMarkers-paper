library(SeuratDisk)
library(Seurat)
library(ggplot2)
setwd("C:/Users/Owner/OneDrive - Johns Hopkins/Desktop/spQSP-HCC/")

#Please note, the single cell data contains samples for all patients
#I separated, and preprocessed all data in scanpy (more-convenient)
#and filtered by nfeature (100), ncell(2), mt_percent(25%)

expression_matrix <- ReadMtx(
  mtx = "singleCell/Cabo_Nivo_aggr_GE/data/hcc13/matrix.mtx", features = "singleCell/Cabo_Nivo_aggr_GE/data/hcc13/genes.tsv",
  cells = "singleCell/Cabo_Nivo_aggr_GE/data/hcc13/barcodes.tsv"
)

sc_hcc13 <- CreateSeuratObject(counts = expression_matrix)
sc_hcc13[['percent_mt']] <- PercentageFeatureSet(sc_hcc13, pattern = '^MT-')
sc_hcc13 <- SCTransform(sc_hcc13, assay = "RNA", vars.to.regress = 'percent_mt',verbose = F)


sc_hcc13 <- GroupCorrelation(sc_hcc13, group.assay = 'RNA', assay = 'SCT', slot = 'scale.data', do.plot = T)


sc_hcc13 <- FindVariableFeatures(sc_hcc13, selection.method = 'vst', nfeatures = 2000)
top10_varfeat <- head(VariableFeatures(sc_hcc13), 10)
top10_varfeat


all_genes <- rownames(sc_hcc13)
sc_hcc13 <- ScaleData(sc_hcc13, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
sc_hcc13 <- RunPCA(sc_hcc13, assay = "SCT", verbose = F)
print(sc_hcc13[['pca']])
DimPlot(sc_hcc13, reduction = 'pca')
#DimHeatmap(sc_hcc13, dims = 1:10)

# determine the dimensionaluty of the dataset
ElbowPlot(sc_hcc13)
# elbow is around PC7-10 -> for HCC02, 10 will capture the true signal

# cluster spots
sc_hcc13 <- FindNeighbors(sc_hcc13, reduction = "pca", dims = 1:10)
sc_hcc13 <- FindClusters(sc_hcc13, verbose = F)

# run UMAP for dimensional reduction

sc_hcc13 <- RunUMAP(sc_hcc13, reduction = "pca", dims = 1:10)
p.umap <- DimPlot(sc_hcc13, reduction = "umap", label = T)


p.umap

### find all markers distinguishing one cluster from all other spots
hcc13_markers <- FindAllMarkers(sc_hcc13, only.pos = T, test.use = "negbinom",
                                min.pct = 0.25, logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc13_markers %>% group_by(cluster) %>% slice_max(n =20, order_by = avg_log2FC))

#write.csv(markers_df,
          #"singleCell/Cabo_Nivo_aggr_GE/data/hcc13/schcc13_markers_cleaned.csv", 
          #row.names = FALSE)
####################################################
#cluster0: cancer cell(ADH1b)
##Merker: "ADH1B"    "CES1"     "APOB"     "NEAT1"    "APOH"  
#         "TF"       "CYP2E1"   "GLUL"     "APOE"     "EPHX1"

#
#cluster1: fibroblasts (with antigen presenting function)
#"TMSB4X"   "HLA-DRA"  "ACTB"     "TMSB10"   "B2M"  
#"IGFBP7"   "CST3"     "HLA-B"    "HSPA1A"   "PTMA"  
#
#cluster2: cancer cell(APOA2)
#"APOA2"    "SERPINA1" "APOC2"    "APOC3"    "APOC1"
#"NME2"     "UQCRQ"    "RBP4"     "HINT1"    "FTH1"
#
####################################################
current.cluster.ids <- c(0, 1, 2)
new.cluster.ids <- c("Cancer cell (ADH1B)", "Immune", "Cancer cell (APOA2)")
sc_hcc13@meta.data$seurat_clusters <- plyr::mapvalues(x = sc_hcc13@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#png(filename = "C:/Users/Owner/OneDrive - Johns Hopkins/Desktop/spQSP-HCC/singleCell/Cabo_Nivo_aggr_GE/figures/heatmap.png", units="in", width=10, height=10, res=600)
DoHeatmap(object = sc_hcc13, 
          features = c("ADH1B",    "CES1",     "APOB",    "NEAT1",    "APOH",  
                       "TF",       "CYP2E1",   "GLUL",     "APOE",     "EPHX1",
                       "TMSB4X",   "HLA-DRA",  "ACTB",     "TMSB10",   "B2M",  
                       "IGFBP7",   "CST3",     "HLA-B",    "HSPA1A",   "PTMA",
                       "APOA2",    "SERPINA1", "APOC2",    "APOC3",    "APOC1",
                       "NME2",     "UQCRQ",    "RBP4",     "HINT1",    "FTH1"),
          group.by = 'seurat_clusters', size=4, angle = 0, hjust = 0.5)
#dev.off()

#png(filename = "C:/Users/Owner/OneDrive - Johns Hopkins/Desktop/spQSP-HCC/singleCell/Cabo_Nivo_aggr_GE/figures/clusters.png", units="in", width=12, height=12, res=600)

DimPlot(sc_hcc13, reduction = "umap",  group.by = 'seurat_clusters') +  ggtitle('')

#dev.off()
###################################################

scale_cogap <- read.csv(file = "singleCell/Cabo_Nivo_aggr_GE/data/cogap/hcc13/projection_15_128.csv", row.names = 1)
scale_cogap <- scale_cogap[rownames(sc_hcc13@meta.data), ]

sc_hcc13@meta.data <- merge(sc_hcc13@meta.data, scale_cogap,
                            by = 'row.names', all = TRUE)

rownames(sc_hcc13@meta.data) <- sc_hcc13@meta.data$Row.names

#png(filename = "C:/Users/Owner/OneDrive - Johns Hopkins/Desktop/spQSP-HCC/singleCell/Cabo_Nivo_aggr_GE/figures/3_pattern_umap.png", units="in", width=12, height=12, res=600)
FeaturePlot(sc_hcc13, features = colnames(scale_cogap)[1:3], ncol = 2) & scale_colour_viridis_c(limits = c(-0.6,0.6), breaks = c(-0.6, 0.0, 0.6))
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(-0.6,0.6), breaks = c(-0.6, 0.0, 0.6)) & scale_fill_viridis_c()

#dev.off()

#png(filename = "C:/Users/Owner/OneDrive - Johns Hopkins/Desktop/spQSP-HCC/singleCell/Cabo_Nivo_aggr_GE/figures/3_Pattern_Cluster.png", units="in", width=8, height=8, res=600)
DimPlot(sc_hcc13, reduction = "umap",  group.by = 'Pattern_Cluster') + ggtitle('Pattern Cluster')
#dev.off()

######################### Atul spaceMarker #########################
p1p8 <- read.table(file='VisiumData/CoGAPS/1541/1541_15Patterns_SpaceMarkers_residual_Pattern_1 x Pattern_8Fri_Mar_18_20.42.04_2022.txt')
p2p8 <- read.table(file='VisiumData/CoGAPS/1541/1541_15Patterns_SpaceMarkers_residual_Pattern_2 x Pattern_8Fri_Mar_18_20.42.04_2022.txt')

DoHeatmap(object = sc_hcc13, 
          features = p2p8$V1,
          group.by = 'Pattern_Cluster', size=4, angle = 0, hjust = 0.5)

############################





