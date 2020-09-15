# Note: This commands following by stLearn_QC.py in QC folder
# You can replace adata by adata_A, adata_B, adata_C or adata_D

st.pp.normalize_total(adata,target_sum=1e4)
st.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
st.pp.scale(adata,max_value=10)

# Create folder tiling before run this command

st.pp.tiling(adata,out_path="../tiling",crop_size =35)
st.pp.extract_feature(adata)	

st.em.run_pca(adata,n_comps=25,random_state=0)
st.spatial.morphology.adjust(adata,use_adata="X_pca",radius=35,method="mean")

st.pp.neighbors(adata,n_neighbors=25,use_rep='X_pca_morphology',random_state=0)
sc.tl.leiden(adata,resolution=0.4,random_state=0,)

st.pl.cluster_plot(adata,use_label="leiden",tissue_alpha=1,spot_size=7,show_legend=True,adata_alpha=1)

sc.tl.umap(adata,min_dist=0.5, spread=3, n_components=2, maxiter=None, alpha=5.0, 
	gamma=5.0, negative_sample_rate=0.5, init_pos='spectral',random_state=0)


sc.pl.umap(adata, color='leiden', add_outline=False, legend_loc='on adata',
           legend_fontsize=0, legend_fontoutline=0,frameon=False,
           title='', wspace=0.4)

