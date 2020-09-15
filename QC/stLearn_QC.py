import stlearn as st
import scanpy as sc

PATH = [Kidney_A_Path, Kidney_B_Path, Kidney_C_Path, Kidney_D_Path]

#########################
#		Kidney A        #
#########################
adata_A = st.Read10X(path=PATH[0])
st.pp.filter_genes(adata_A,min_cells=3)

adata_A.var['mt'] = adata_A.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'

sc.pp.calculate_qc_metrics(adata_A, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata_A, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata_A, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_A, x='total_counts', y='n_genes_by_counts')

adata_A = adata_A[adata_A.obs.n_genes_by_counts < 4500, :]

#########################
#		Kidney B        #
#########################

adata_B = st.Read10X(path=PATH[1])
st.pp.filter_genes(adata_B,min_cells=3)

adata_B.var['mt'] = adata_B.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata_B = adata_B[:, ~adata_B.var['mt']]

sc.pp.calculate_qc_metrics(adata_B, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata_B, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata_B, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_B, x='total_counts', y='n_genes_by_counts')

adata_B = adata_B[adata_B.obs.n_genes_by_counts < 3500, :]

#########################
#		Kidney C        #
#########################

adata_C = st.Read10X(path=PATH[2])
st.pp.filter_genes(adata_C,min_cells=3)
adata_C.var['mt'] = adata_C.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata_C = adata_C[:, ~adata_C.var['mt']]

sc.pp.calculate_qc_metrics(adata_C, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata_C, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata_C, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_C, x='total_counts', y='n_genes_by_counts')

adata_C = adata_C[adata_C.obs.n_genes_by_counts < 3000, :]

#########################
#		Kidney D        #
#########################

adata_D = st.Read10X(path=PATH[3])
st.pp.filter_genes(adata_D,min_cells=3)
adata_D.var['mt'] = adata_D.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'

st.add.add_loupe_clusters(adata_D,loupe_path=DEMULTIPLEX_BARCODES_PATH,key_add="Tissue")
adata_D = adata_D[adata_D.obs["Tissue"]=="Left"]
sc.pp.calculate_qc_metrics(adata_D, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata_D, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata_D, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_D, x='total_counts', y='n_genes_by_counts')