# Note: This commands following by SPOTlight_deconvolution.R in Deconvolution folder
# You can replace adata by adata_A, adata_B, adata_C or adata_D

import stlearn as st

adata = st.Read10X(path=ST_SEQ_PATH)
st.add.auto_annotate(adata,annotation_path="Kidney_Lake_deconvolution_result.csv")

st.pl.deconvolution_plot(adata,cmap="tab10",threshold=0.0,spot_size=12)