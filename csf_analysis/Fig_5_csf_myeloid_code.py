# Python script showing the Myeloid CSF cNMF analysis
# Regan Murphy 

#Imports used below

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sklearn
import scipy
import h5py
import anndata as ad

# for cosine similarity
from sklearn.metrics.pairwise import cosine_similarity

import re
from sklearn.linear_model import LinearRegression
from scipy import stats
from itertools import combinations

'''
This file walks through the process of preparing the single cell object for cNMF and then downstream analysis of the output matrices
'''

#obj_filtered_atlas = anndata object containing myeloid cells from the CSF
#                     equivalent to a seurat object
#                     should have raw counts saved as object.X
#                     already QC filtered
#                     can contain non-myeloid cells

# Step 1 - Identification of Myeloid Cells

# 1.1 Normalization
# Saving count data
obj_filtered_atlas.layers["counts"] = obj_filtered_atlas.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(obj_filtered_atlas)

# Logarithmize the data
sc.pp.log1p(obj_filtered_atlas)

# 1.2 Selecting Variable Genes for Leiden analysis
obj_filtered_atlas.obs.index = obj_filtered_atlas.obs.index.astype(str)

sc.pp.highly_variable_genes(obj_filtered_atlas, n_top_genes=4000, batch_key="sample")
sc.pl.highly_variable_genes(obj_filtered_atlas)

# 1.3 Dimensionality Reduction
sc.tl.pca(obj_filtered_atlas)
sc.pl.pca_variance_ratio(obj_filtered_atlas, n_pcs=50, log=True)

# 1.4 UMAP
sc.pp.neighbors(obj_filtered_atlas, n_pcs=20) #select n_pcs based on the variance_ratio plot
sc.tl.umap(obj_filtered_atlas)

# 1.5 Clustering & Identifying Myeloid Clusters

## As needed, adjust the resolution of clusters (higher = more granular), recommended to run a few resolutions
sc.tl.leiden(obj_filtered_atlas, flavor="igraph", n_iterations=2, resolution = 0.1)

# Markers used for annotation
# IF errors thrown about genes not present in matrix, remove from this dictionary & re-run
corrected_myeloid_markers = {
    "immune": ["PTPRC"],
    "oligo": ["THY1", "PLP1", "APOD"],
    "pericyte": ["PDGFRB", "MCAM"],
    "endothelial": ["PECAM1", "ACTA2"],
    "cDC": ["HLA-DRA", "HLA-DRB1","HLA-DPB1", "AREG", "FCER1A","LAMP3", ],
    "neutrophil": ["S100A8", "IFITM2", "FCGR3B"],
    "monocyte": ["LYZ", "VCAN", "FN1"],
    "mast": [ "KIT"], #"MS4A2",
    "myeloid": ["MPO","ITGAM","ITGAX","GPNMB", "TMEM119", "MRC1", "CD163", "IL1B", "CD83", "CCL3", "CD68","FUT4"],
    "t_cell": ["CD3E"],
    "b_cell": ["MS4A1"],
    "plasma": ["SDC1"],
    "cancer": ["EGFR","PTPRZ1","CSRP2","SOX2"]#"AQP4"
}

sc.pl.dotplot(obj_filtered_atlas, corrected_myeloid_markers, groupby='leiden')

# Identify which clusters show myeloid markers & subset the obj_filtered_atlas
barcodes_myeloid_rough = obj_filtered_atlas.obs[
    obj_filtered_atlas.obs["leiden"].isin(['0', '1','2','3','4','5','6','7','8','9','10','11','12','13','14'])
].index.tolist()
len(barcodes_myeloid_rough)

obj_subset_myeloid_rough = obj_filtered_atlas[barcodes_myeloid_rough]

obj_subset_myeloid_rough

# If clustering isn't completely clean, reclustering the initial selection would follow the same steps

barcodes_myeloid = obj_subset_myeloid_rough.obs[
    obj_subset_myeloid_rough.obs["leiden"].isin(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
 "11", "12", "13", "14", "15", "16", "17", "18", "19",
 "20", "21", "23"])
].index.tolist()
len(barcodes_myeloid)

obj_subset_myeloid = obj_subset_myeloid_rough[barcodes_myeloid]

obj_subset_myeloid


# selection of 2000 variable genes
# Function used to select variable genes
def select_highly_variable_genes(adata_raw_counts, norm_target_sum, f1_mean_exp_cutoff, f2_exp_threshold, f2_fraction_total_cells, num_var_genes):
    '''
    This function performs a highly variable gene selection with two preliminary filter criteria: 
        (1) mean expression cutoff, (2) per-cell expression threshold in specified fraction of the total dataset. 
        Of these filtered genes, a specified number of highest variable genes is picked and the adata object is subset accordingly. 
    
    Args: 
        adata_raw_counts: an adata object with only raw counts (before normalization)
        norm_target_sum: the normalization target sum (default: 1e4)
        f1_mean_exp_cutoff: the mean expression threshold for filter 1 (default: 0.01)
        f2_exp_threshold: the per-cell expression threshold for filter 2 (default: 0.1)
        f2_fraction_total_cells: the fraction of total cells to have the f2_exp_threshold (default: 0.001)
        num_var_genes: the number of highly variable genes to pick (default: 2000)
        
    Returns:
        an adata object (raw counts) with only the top variable genes
    '''
    
    
    ######################
    # Normalize Raw Data #
    ######################
    
    adata_normalized = adata_raw_counts.copy()
    
    # Normalizing to median total counts
    sc.pp.normalize_total(adata_normalized, 
                          target_sum = norm_target_sum)

    # Logarithmize the data
    sc.pp.log1p(adata_normalized)

    
    ####################################
    # Filter 1: Mean Expression Cutoff #
    ####################################
    
    # Compute mean expression for each gene across all cells
    mean_expression = np.mean(adata_normalized.X.toarray(), axis=0)  # Convert sparse matrix to dense

    # Get the indices of genes with mean expression > f1_mean_exp_cutoff
    f1_genes_to_keep = mean_expression > f1_mean_exp_cutoff

    # Extract gene names from the var DataFrame
    f1_selected_genes = adata_normalized.var_names[f1_genes_to_keep]

    # Display the selected genes
    print("Filter 1")
    print(f"Number of genes with mean expression > 0.01: {len(list(f1_selected_genes))}")
    
    #############################################################
    # Filter 2: High Expression in Pct of Total Cells Threshold #
    #############################################################
    
    threshold_cells = int(f2_fraction_total_cells * adata_normalized.n_obs)

    # Identify number of cells with expression of genes greater than "f2_exp_threshold"  
    cell_counts = np.sum(adata_normalized.X > f2_exp_threshold, axis=0).A1  # Count cells where expression > 0.1

    # Identify genes expressed in at least `threshold_cells`
    f2_genes_to_keep = cell_counts >= threshold_cells

    # Extract gene names that meet the condition
    f2_selected_genes = adata_normalized.var_names[f2_genes_to_keep]

    # Display the selected genes
    print("\nFilter 2")
    print(f"Number of genes with > {f2_exp_threshold} expression in at least {threshold_cells} cells : {len(list(f2_selected_genes))}")
    
    
    #############################################
    # Identify the genes that meet both filters #
    #############################################
    
    filtered_genes = list(set(f1_selected_genes).intersection(set(f2_selected_genes)))

    print(f"\nNumber of genes meeting both filters : {len(list(filtered_genes))}")

    
    ################################
    # Pick the top 'num_var_genes' #
    ################################
    
    # subset the raw counts
    adata_raw_counts_filtered = adata_raw_counts[:, filtered_genes] 
    
    # identify the top 'num_var_genes' variable genes
    sc.pp.highly_variable_genes(adata_raw_counts_filtered, 
                                   flavor = "seurat_v3", 
                                   n_top_genes = num_var_genes)
    
    df_highly_var = (pd.DataFrame(adata_raw_counts_filtered.var["highly_variable"]))

    top_var_genes = list(df_highly_var[df_highly_var["highly_variable"] == True].index)
    
    # subset the adata object to only the top variable genes
    return adata_raw_counts_filtered[:, top_var_genes]

# read in your full single cell object (raw counts, qc filtered)
# subset to barcodes_myeloid as shown aboveto get unprocessed_obj_myeloid_subset

unprocessed_obj_myeloid_subset_filtered = select_highly_variable_genes(adata_raw_counts = unprocessed_obj_myeloid_subset, 
                                                                                             norm_target_sum = 1e4, 
                                                                                             f1_mean_exp_cutoff = 0.01, 
                                                                                             f2_exp_threshold = 0.1, 
                                                                                             f2_fraction_total_cells = 0.001, 
                                                                                             num_var_genes = 2000 )

unprocessed_obj_myeloid_subset_filtered

# you should end with an object that has ~80-90k myeloid cells (depending on how easy cluster-selection was)
#                                       and 2000 genes
#                                       Important -- still needs raw counts 

# cNMF running & program annotation (terra workflows)
# Dylan Kotliar's cnmf script url: https://github.com/dylkot/cNMF
'''
Params: adata_file = saved unprocessed_obj_myeloid_subset_filtered
        output_directory = where you want it to save outputs
        density_threshold = 0.02
        extra_disc_space = 75
        k_range = ["4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26" ]
        memory = 64G
        num_genes = 2000
        num_iter = 500
        seed = 42
        selected_k 
    
        First run of the workflow should leave selected_k blank, it will output an error and stability plot across the k-range
        Our first run we picked k = 14 for a local maximum of the stability 
'''
# Ty Miller Lab cnmf_post_processing_program_annotation script url: https://github.com/cwru-non-academic/bioinformatics_workflows/blob/d70c679dce9b77cd552022d451b22fafb39c7950/scripts/cnmf_gprofiler_excel_local_paths.py
# This workflow is essentially a wrapper for the gProfiler API to do GSEA using the top 100 genes
# per program as ranked by spectra
'''
Params: gene_count = 100
        output_file_name = file name
        remote_gene_spectra_score_mtx_path = /../cnmf_out_dir/cnmf_run.gene_spectra_score.k_{selected_k}.dt_0_02.txt
        remote_gmt_path = location of GMT annotation file
        remote_output_directory_path = where it will be saved
'''

# After running the annotation script you get a spreadsheet with the matches in the gmt file
# ranked by the p-values with the genes ranked by spectra score for each program
# highlighted are any matches in the gmt file generated at MGH (by Tyler Miller et al, Nature 2025)

# If your myeloid selection was clean there should be only myeloid based annotations
# If there were some weird cells there may be some programs that show other cell types
# Those cells can be filtered out by usage (look below for generating the matrices) and then cnmf can be reran with cleaner cells
# We ran cNMF twice k = 14 & k = 13

# creation of usage matrix & variations (explain differences)

# From the cNMF output folder look for the below file name and load it in as a pandas DataFrame
# df_usage_matrix_k13 = cNMF_output_dir/cnmf_run.usages.k_{selected_k}.dt_0_02.consensus.txt

# load in the sample metadata as a pandas DataFrame
#df_metadata = '/home/resources/...'

df_metadata.index.name = "sample"
df_metadata

# insert the program names that were identified in the gProfiler excel analysis

# EXAMPLE BELOW
prog_cols_k13 = ["cDC1",
"Suppressive Monocyte",
"CD14+ Monocyte",
"cDC3",
"pDC",
"p6-tcell",
"Interferon",
"Scavenger Macrophage",
"CSF Inflammatory",
"Inflammatory Macrophage",
"Mixed Myeloid",
"Complement",
"cDC2"
]

barcode_len = 17
usage_threshold = 20


colnames = ['Barcodes'] + prog_cols_k13
df_usage_matrix_k13.columns = colnames
df_usage_matrix_k13 = df_usage_matrix_k13.set_index('Barcodes')
df_usage_matrix_k13.head()

df_usage_matrix_k13 = df_usage_matrix_k13.drop('p6-tcell', axis=1) #this program had low usage and was dropped after cNMF
df_usage_matrix_k13

# myeloid_h5ad = load in the single-cell object 

# merges in any cell level metadata (like sample label) based on the barcodes
df_usage_matrix_k13 = df_usage_matrix_k13.merge(
    myeloid_h5ad.obs[['sample']],          
    left_on=df_usage_matrix_k13.columns[0],     
    right_index=True,                
    how='left'                        
)
# from the sample label generate patient labels
df_usage_matrix_k13['patient'] = df_usage_matrix_k13['sample'].str.split('D', n=1).str[0]

# merge in sample level metadata
df_usage_matrix_k13 = df_usage_matrix_k13.merge(
    df_metadata,
    left_on='patient',
    right_index=True,
    how='left'
)

# extract timepoint labels from sample
df_usage_matrix_k13['Day'] = df_usage_matrix_k13['sample'].str.extract(r'(D\d+)')[0]

# update the prog_cols_k13 (program column names) after dropping p6-tcell
prog_cols_k13 = ["cDC1",
"Suppressive Monocyte",
"CD14+ Monocyte",
"cDC3",
"pDC",
"Interferon",
"Scavenger Macrophage",
"CSF Inflammatory",
"Inflammatory Macrophage",
"Mixed Myeloid",
"Complement",
"cDC2"]

# row_normalize (ie. the usage of programs across each cell will add to 100%)
data_cols = df_usage_matrix_k13[prog_cols_k13]
norm_sums = data_cols.sum(axis=1).tolist()
df_per_cell_usage_matrix_k13 = data_cols.div(norm_sums, axis=0).mul(100)


# Annotate cells with their most used identity program and most used activity program
id_prog = ["cDC1",
"Suppressive Monocyte",
"CD14+ Monocyte",
"cDC3",
"pDC",
"Inflammatory Macrophage",
"Mixed Myeloid",
"cDC2"]
act_prog = ["Interferon",
"Scavenger Macrophage",
"CSF Inflammatory",
"Complement"]

backup_id = 'activity'
backup_act = 'identity'
id_max_vals = df_per_cell_usage_matrix_k13[id_prog].max(axis=1)
id_max_idx = df_per_cell_usage_matrix_k13[id_prog].idxmax(axis=1)
id_primary_program = id_max_idx.where(id_max_vals >= usage_threshold, backup_id)
df_per_cell_usage_matrix_k13['Primary Identity'] = id_primary_program

act_max_vals = df_per_cell_usage_matrix_k13[act_prog].max(axis=1)
act_max_idx = df_per_cell_usage_matrix_k13[act_prog].idxmax(axis=1)
act_primary_program = act_max_idx.where(act_max_vals >= usage_threshold, backup_act)
df_per_cell_usage_matrix_k13['Primary Act'] = act_primary_program

pd.DataFrame(df_per_cell_usage_matrix_k13[["Primary Identity", "Primary Act"]].value_counts()).head(60)

# Generate the per-sample usage matrix
df_cell_counts_per_sample_k13 = pd.DataFrame(df_usage_matrix_k13["sample"].value_counts())

# collapse down to the sum of cells with > 20% usage of each program on a per-sample basis
df_sample_usage_matrix_k13 = (df_per_cell_usage_matrix_k13[prog_cols_k13].gt(usage_threshold).astype(int).join(df_usage_matrix_k13['sample']).groupby('sample').sum())


patient_prog_sums_k13 = df_cell_counts_per_sample_k13["count"]
df_patient_norm_usage_matrix_k13 = df_sample_usage_matrix_k13.div(patient_prog_sums_k13, axis=0).mul(100)

# Append sample-level metadata
df_sample_metadata = df_usage_matrix_k13.loc[:, ~df_usage_matrix_k13.columns.isin(prog_cols_k13)]
df_sample_metadata = df_sample_metadata.groupby("sample", as_index=True).first()

df_patient_norm_usage_matrix_k13_with_metadata = df_patient_norm_usage_matrix_k13.copy()
df_patient_norm_usage_matrix_k13_with_metadata = df_patient_norm_usage_matrix_k13_with_metadata.merge(df_sample_metadata, on='sample')

# Subsetting to high-dose patients only
high_dose_pts = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12']
low_dose_pts = ['P13', 'P14', 'P15', 'P16', 'P17', 'P18']

high_dose_usage = df_usage_matrix_k13[df_usage_matrix_k13['patient'].isin(high_dose_pts)]


# row_normalize (ie. the usage of programs across each cell will add to 100%)
hd_data_cols = high_dose_usage[prog_cols_k13]
hd_norm_sums = hd_data_cols.sum(axis=1).tolist()
hd_df_per_cell_usage_matrix_k13 = hd_data_cols.div(hd_norm_sums, axis=0).mul(100)

# Quality checking the matrix if there are any divide by 0 --> NaN errors
nan_barcodes = hd_df_per_cell_usage_matrix_k13[hd_df_per_cell_usage_matrix_k13.isna().any(axis=1)].index
hd_df_per_cell_usage_matrix_k13 = hd_df_per_cell_usage_matrix_k13.drop(index=nan_barcodes)

hd_df_cell_counts_per_sample_k13 = pd.DataFrame(high_dose_usage["sample"].value_counts())

# collapse down to the sum of cells with > 20% usage of each program on a per-sample basis
hd_df_sample_usage_matrix_k13 = (hd_df_per_cell_usage_matrix_k13[prog_cols_k13].gt(usage_threshold).astype(int).join(high_dose_usage['sample']).groupby('sample').sum())

hd_patient_prog_sums_k13 = hd_df_cell_counts_per_sample_k13["count"]
hd_df_patient_norm_usage_matrix_k13 = hd_df_sample_usage_matrix_k13.div(hd_patient_prog_sums_k13, axis=0).mul(100)
hd_df_patient_norm_usage_matrix_k13 = hd_df_patient_norm_usage_matrix_k13.dropna()

hd_df_patient_norm_usage_matrix_k13_with_metadata = hd_df_patient_norm_usage_matrix_k13.copy()
hd_df_patient_norm_usage_matrix_k13_with_metadata = hd_df_patient_norm_usage_matrix_k13_with_metadata.merge(hd_df_sample_metadata, on='sample')

# Creating the per-cell usage heatmap 

# Assigning an interpretable order for the programs, grouping by activity and identity
ordered_prog_cols_k13 = ["CSF Inflammatory","Complement","Interferon",
"Scavenger Macrophage","Inflammatory Macrophage",
"Suppressive Monocyte","CD14+ Monocyte",
"Mixed Myeloid","cDC1","cDC2", "cDC3","pDC"]

hd_ordered_df = hd_df_per_cell_usage_matrix_k13.loc[:, ordered_prog_cols_k13]
hd_ordered_df_T = hd_ordered_df.T
hd_ordered_df_T

meta_df = high_dose_usage[['patient', 'Day', 'Reponder']]

meta_df=meta_df.loc[hd_ordered_df_T.columns]

patient_palette = sns.color_palette("tab20", meta_df['patient'].nunique())
patient_lut = dict(zip(meta_df['patient'].unique(), patient_palette))
patient_colors = meta_df['patient'].map(patient_lut)

day_lut = {
    'D0': '#6c0cf2',
    'D7': '#f26b0c',
    'D21': '#0cf26c'
}
day_colors = meta_df['Day'].map(day_lut)

responder_lut = {
    True: '#1b9e77',
    False: '#d95f02'
}
responder_colors = meta_df['Reponder'].map(responder_lut)

col_colors = pd.DataFrame({
    'Patient': patient_colors,
    'Day': day_colors,
    'Responder': responder_colors
}, index=meta_df.index)

g = sns.clustermap(
    hd_ordered_df_T,
    cmap="viridis",
    row_cluster=False,
    col_cluster=True,
    figsize=(12, 8),
    xticklabels=False,
    yticklabels=True,
    col_colors=col_colors
)

# DataFrame manipulation to prep for other downstream analyses

# long_df: for each sample (unique patient and timepoint) have a row for every program with its usage

hd_timepoint_df = hd_df_patient_norm_usage_matrix_k13_with_metadata.copy()
hd_long_df = hd_timepoint_df.melt(
    id_vars=['patient', 'Day', 'Subject', 'Age', 'Sex', 'Dose Level', 'Dose (cells)', 'DoseClass',
       'Reponder', 'PFS_Final_2025_11_14', 'OS_Final_2025_11_14','OS_Final_2025_11_12 FROM ZEV',
       'max tumor regression %', 'cohort'],# identifiers to keep
    var_name='Program',                    # new column name for program names
    value_name='Usage'                     # new column name for usage values
)

hd_df_long = hd_long_df.sort_values(['Program', 'patient', 'Day'])

# wide_df: for every patient create a column for each program's usage at each time point
#           additionally calculate all the program usage time point differences

df_wide = hd_df_patient_norm_usage_matrix_k13_with_metadata.pivot(
    index="patient",
    columns="Day",
    values=prog_cols_k13
)

df_wide.columns = [f"{prog}_{day}" for prog, day in df_wide.columns]
df_wide = df_wide.reset_index()

patient_meta = hd_df_patient_norm_usage_matrix_k13_with_metadata.drop_duplicates("patient")[
    ["patient", "Subject", "Age", "Sex", "Dose Level", 
     "Dose (cells)", "DoseClass", "Reponder",
     "PFS_Final_2025_11_14", "OS_Final_2025_11_14",
     "max tumor regression %", 'OS_Final_2025_11_12 FROM ZEV', 'cohort']
]
df_wide = df_wide.merge(patient_meta, on="patient", how="left")

time_pairs = [
    ("D0", "D7"),
    ("D7", "D21"),
    ("D0", "D21")
]

for t1, t2 in time_pairs:
    for prog in prog_cols_k13:
        col1 = f"{prog}_{t1}"
        col2 = f"{prog}_{t2}"
        df_wide[f"{prog}_{t2}_minus_{t1}"] = df_wide[col2] - df_wide[col1]
