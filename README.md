# SpatialKidney

## stLearn analysis pipeline:

- Step 1 - QC for raw data: ./QC/stLearn_QC.py
- Step 2 - Find cell types at cluster level: ./CellTypeAnnotation/stlearn_Clustering.py
- Step 3 - Compare gene expression: ./DifferentialExpression/stLearn_DE.py
- Step 4 - Find cell types at spot level: ./Deconvolution/SPOTlight_deconvolution.R and Deconvolution/visualize_deconvolution.py
- Step 5 - Predict cell-to-cell interactions: ./Cell_Cell_Interactions/
- Step 6 - Investigate kidney disease associated SNPs and genes: ./SpatialGenomics/SNP_Analysis.sh  


## Datasets:

- Link 1 human kidney  -> https://cloud.rdm.uq.edu.au/index.php/s/rqnfW4a5NdxqyKW -> password: contact us to get the password 
- Link 2  human kidney  -> https://cloud.rdm.uq.edu.au/index.php/s/DxPY9sxTR27yXLQ -> password: contact us to get the password 
- Link 3  human kidney  -> https://cloud.rdm.uq.edu.au/index.php/s/m3Dgj29QetZpLmJ -> password: contact us to get the password 
- Link 4  Legacy LP3 (mouse adult and human kidney -> https://cloudstor.aarnet.edu.au/plus/s/NkrQvVjsTtWvKiU) -> no password
