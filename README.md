# SpatialKidney

## stLearn analysis pipeline:

- Step 1 - QC for raw data: ./QC/stLearn_QC.py
- Step 2 - Find cell types at cluster level: ./CellTypeAnnotation/stlearn_Clustering.py
- Step 3 - Compare gene expression: ./DifferentialExpression/stLearn_DE.py
- Step 4 - Find cell types at spot level: ./Deconvolution/SPOTlight_deconvolution.R and Deconvolution/visualize_deconvolution.py
- Step 5 - Predict cell-to-cell interactions: ./Cell_Cell_Interactions/
- Step 6 - Investigate kidney disease associated SNPs and genes: ./SpatialGenomics/SNP_Analysis.sh  
