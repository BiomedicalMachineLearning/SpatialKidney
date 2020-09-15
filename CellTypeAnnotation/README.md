# stLearn and Seurat strategies for cell type annotation

## SeuratSpatialAnalysisPipeline.Rmd

This R notebook will analyse a sample from SpaceRanger output (containing `spatial` directory and `filtered_feature_bc_matrix.h5` file), through QC, dimensionality reduction, clustering and marker detection. The user must edit the first code chunk to provide the path to this SpaceRanger output and a directory to store the output files. This R notebook can be converted to an R script by calling:

```
knitr::purl("/path/to/SeuratSpatialAnalysisPipeline.Rmd")
# will save a file /path/to/SeuratSpatialAnalysisPipeline.R
```
