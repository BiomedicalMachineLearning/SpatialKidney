# stLearn and Seurat strategies for cell type annotation

## SeuratSpatialAnalysisPipeline.Rmd

This R notebook will analyse a sample from SpaceRanger output (containing `spatial` directory and `filtered_feature_bc_matrix.h5` file), through QC, dimensionality reduction, clustering and marker detection. The user must edit the first code chunk to provide the path to this SpaceRanger output and a directory to store the output files. This R notebook can be converted to an R script by calling:

```
knitr::purl("/path/to/SeuratSpatialAnalysisPipeline.Rmd")
# will save a file /path/to/SeuratSpatialAnalysisPipeline.R
```

Users wishing to run the cell cycle prediction step will need to download the `mouse.cc.genes.Rdata` file and update the path in the `load("/path/to/mouse.cc.genes.Rdata")` line. 

## SeuratSpatialAnalysisPipeline.Rmd

This R notebook will perform label transfer to annotate a query dataset with annotations from a reference dataset. The user must edit the first code chunk to provide the path to the query and reference RDS files (must be Seurat objects). This R notebook can be converted to an R script by calling:

```
knitr::purl("/path/to/SeuratSpatialAnalysisPipeline.Rmd")
# will save a file /path/to/SeuratSpatialAnalysisPipeline.R
```
