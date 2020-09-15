Information about reference datasets used for label transfer

# Liao et al

Liao, J., Yu, Z., Chen, Y. et al. Single-cell RNA sequencing of human kidney. Sci Data 7, 4 (2020). https://doi.org/10.1038/s41597-019-0351-8

Data availability: GEO database, accession number GSE131685

Code adapted from: https://github.com/lessonskit/Single-cell-RNA-sequencing-of-human-kidney

Files provided here:

* `Liao_AnalyseData.Rmd` - code adapted from original publication. Running the authors' code resulted resulted in 11 clusters, while 10 were presented in the original paper (most likely due to subsequent updates to the Seurat algorithm). The cluster counts in Figure 1C of the original publication (see file `Liao_Fig1CCounts`) were used to convert between the clusters found by Liao et al. and the clusters we found (see file `Liao_ClusterConversion`).
* `Liao_SimplifyAnnotations.Rmd` - annotations for several cell types were combined to provide a simplified annotation labelling the broad structural components of the kidney
* `Liao_Metadata.txt` - metadata file listing all cells present in the dataset, the clusters we identified, the original specific annotations, and the simplified annotations used for label transfer.
