library(Seurat)
library(SPOTlight)

# You can change the name lake or liao depend on your reference dataset. Here is the example of deconvolution using Lake dataset
# The reference dataset could be found in the folder ReferenceDatasets

lake <- readRDS(LAKE_REFERENCE_SEURAT_OBJECT_PATH)
lake[["RNA"]] <- lake[["lake"]]
lake@active.assay <- "RNA"	
lake <- subset(lake, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)
lake <- SCTransform(lake, verbose = FALSE)
lake <- RunPCA(lake, dims=1:30, verbose = FALSE)
lake <- RunUMAP(lake, dims = 1:30, verbose = FALSE)
lake <- FindNeighbors(lake, dims = 1:30, verbose = FALSE)
lake <- FindClusters(lake, verbose = FALSE)


kidney <- Load10X_Spatial(ST_SEQ_PATH)
kidney <- SCTransform(kidney, assay = "Spatial", verbose = FALSE)


Idents(object = lake) <- lake@meta.data$celltype_simple
cluster_markers_all <- Seurat::FindAllMarkers(object = lake, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE, 
                                              , min.pct = 0.5, logfc.threshold = 0.5)


set.seed(123)
spotlight_ls <- spotlight_deconvolution(se_sc = lake,
                                      counts_spatial = kidney@assays$Spatial@counts,
                                      clust_vr = "celltype_simple",
                                      cluster_markers = cluster_markers_all,
                                      cl_n = 50, # 100 by default
                                      hvg = 3000,
                                      ntop = NULL,
                                      transf = "uv",
                                      method = "nsNMF",
                                      min_cont = 0.09)

decon_mtrx <- spotlight_ls[[2]]
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

kidney@meta.data <- cbind(kidney@meta.data, decon_mtrx)

tmp = as.data.frame(decon_mtrx)
row.names(tmp) <- row.names(kidney@meta.data)
write.csv(t(tmp[1:(length(tmp)-1)]),"Kidney_Lake_deconvolution_result.csv")

