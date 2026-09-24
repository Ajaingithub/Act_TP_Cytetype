#region Day1
import anndata
import scanpy as sc
from cytetype import CyteType
import mygene
import pandas as pd
import os

# ------ Example Scanpy Pipeline ------
#  Skip this step if you already have clusters and marker genes in an AnnData object. 
adata = anndata.read_h5ad("/mnt/data/projects/.immune/Mayo/Activation_TP/processed_data/CD4_D1.h5ad")
# Since the log normalization is already done
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata_full)

(day1.obs['clusters'].astype('str') + "-" +  day1.obs['cytetype_annotation_clusters'].astype('str'))

adata_full = adata.raw.to_adata()
### Changing the gene name
mg = mygene.MyGeneInfo()
ensg_ids = adata_full.var_names.tolist()
gene_info = mg.querymany(ensg_ids, scopes='ensembl.gene', fields='symbol', species='human')
gene_map = pd.DataFrame(gene_info)[['query', 'symbol']].dropna()

adata_full.var['gene_symbol'] = adata_full.var_names.map(
    dict(zip(gene_map['query'], gene_map['symbol']))
)

adata_full2 = adata_full[:,~adata_full.var['gene_symbol'].isna()].copy() # remove Nan values
adata_full3 = adata_full2[:,~adata_full2.var['gene_symbol'].duplicated()] # remove duplicated
adata_full3.var.index = adata_full3.var['gene_symbol']

# sc.pp.highly_variable_genes(adata_full3, n_top_genes=2000)
# sc.pp.pca(adata_full3)
# sc.pp.neighbors(adata_full3)
# sc.tl.leiden(adata_full3, key_added="clusters")

adata_full3.obs['clusters'] = adata_full3.obs['seurat_clusters'].astype(str).str.replace("Cluster ","")
sc.tl.rank_genes_groups(adata_full3, groupby="clusters", method="t-test")

adata_full3.var['gene_symbols'] = adata_full3.var.index
annotator = CyteType(adata_full3, group_key="clusters")
adata_full3 = annotator.run(
    study_context="CD4+ T cells extracted from PBMC of healthy individuals stimulated for 24 hours with CD3-CD28 beads",
)

os.chdir(savedir)
sc.pl.umap(adata_full3,color="cytetype_annotation_clusters",size=5, save = "_cytetype_annotation.pdf")
sc.pl.umap(adata_full3,color="clusters",size=5, save = "_clusters.pdf")
sc.pl.umap(adata_full3,color="seurat_clusters",size=5, save = "_seurat_clusters.pdf")
sc.pl.umap(adata_full3,color="samples",size=4, save = "_sample.pdf")

### Lets try to do annotation based on the seurat clusters
sc.tl.rank_genes_groups(adata_full3, groupby="seurat_clusters", method="t-test")
adata_full3.var['gene_symbols'] = adata_full3.var.index
annotator = CyteType(adata_full3, group_key="clusters")
adata_full3 = annotator.run(
    study_context="CD4+ T cells extracted from PBMC of healthy individuals",
)

sc.pl.umap(adata_full3,color="cytetype_annotation_clusters",size=4, save = "_cytetype_annotation_seurat_clusters.pdf")

markers = [
    "LEF1","TCF7","CCR7","BCL2","SELL","BACH2","FOXP3","IKZF1","HSP90AA1","ICOS","GIMAP7","BCL11B","CD40LG","PTPRC",
    "PCNA","MKI67","TOP2A","STAT1","EGR3","GBP5","GBP2","B2M","IL7R","PIK3CD","RPS29","RPL21","RPS8","FYN","IRF4","CD69"
    ]

sc.pl.dotplot(adata_full3,markers,'cytetype_annotation_clusters', standard_scale = "var",save = "_cytetype_annotation.pdf")

markers = [
    "TCF7","CCR7","GIMAP7","IL7R","BCL11B","SELL","CD40LG","HSP90AA1","HSP90AB1",
    "HSPD1","ARHGAP15", "SESN3", "INPP4B","PCNA","MKI67","TOP2A","RPS29",
    "RPL21","RPS8","FYN","IRF4","CD69","CTLA4","PDCD1","TOX",
    ]

sc.pl.dotplot(adata_full3,markers,'cytetype_annotation_clusters', standard_scale = "var",save = "_cytetype_annotation_2.pdf")
sc.pl.heatmap(adata_full3,markers,'cytetype_annotation_clusters',standard_scale = "var",save = "_cytetype_annotation_2.pdf")

datadir = "/mnt/data/projects/.immune/Mayo/Activation_TP/processed_data/"
savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/"
adata_full3.write_h5ad(os.path.join(savedir,"CD4_Day1_annotation.h5ad"))

sc.tl.rank_genes_groups(adata_full3, groupby="cytetype_annotation_clusters", method="t-test", pts=True)

# Extract results
res = adata_full3.uns['rank_genes_groups']
groups = res['names'].dtype.names  # all cell type names

# Define save directory
savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/Table/DE_ttest/"
os.makedirs(savedir, exist_ok=True)

# Loop through all cell types
for group in groups:
    genes = res['names'][group]
    # build DataFrame
    de_df = pd.DataFrame({
        'gene': genes,
        'logfoldchange': res['logfoldchanges'][group],
        'pvals': res['pvals'][group],
        'pvals_adj': res['pvals_adj'][group],
    })
    # add expression percentages (if pts=True was used)
    if 'pts' in res:
        pts = pd.DataFrame(res['pts'], index=adata_full3.var_names)
        pts_rest = pd.DataFrame(res['pts_rest'], index=adata_full3.var_names)  
        # subset to genes in this group
        pts_group = pts.loc[genes]
        pts_rest_group = pts_rest.loc[genes]
        de_df['pct_expr_in_group'] = pts_group[group].values
        de_df['pct_expr_in_rest'] = pts_rest_group[group].values
    # save to CSV
    group=group.replace("/","_")
    out_path = os.path.join(savedir, f"{group}_DE_results.csv")
    de_df.to_csv(out_path, index=False)
    print(f"Saved {group} markers → {out_path}")

#region quantification
import pandas as pd
import os
savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/"

df = adata_full3.obs[['samples','cytetype_annotation_clusters']].copy()
counts = df.groupby(['samples','cytetype_annotation_clusters']).size().reset_index(name='cell_count')
counts['Age'] = counts['samples'].astype(str).str.replace(r'(\D+)\d.*', r'\1', regex=True)

counts['total_cells'] = counts.groupby('samples')['cell_count'].transform('sum')
counts['percent'] = (counts['cell_count'] / counts['total_cells']) * 100

counts['Age'] = pd.Categorical(counts['Age'],categories=['HY','HO','FO'],ordered=True)

## Just checking
counts.loc[counts['samples'].isin(["HY4D1N4"]),"percent"].sum()

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 8))
sns.barplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    palette='Set2',
    errorbar='sd'
)

sns.stripplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    dodge=True,
    palette='Set2',
    linewidth=0.8,
    edgecolor='black',
    size=5,
    marker='o',
    alpha=0.7,
    legend=False  # keep legend from barplot only
)

### This is T-test if data is normally distributed
# from scipy.stats import ttest_ind
# pvals = {}
# for ct, subdf in counts.groupby('identity'):
#     grp = subdf.groupby('response')['percent']
#     if len(grp) == 2:
#         r1, r2 = [v.values for k, v in grp]
#         # two-sided independent t-test
#         stat, p = ttest_ind(r1, r2, equal_var=False)  # Welch’s t-test (safer)
#         pvals[ct] = p

from scipy.stats import f_oneway
pvals = {}
for ct, subdf in counts.groupby('cytetype_annotation_clusters'):
    groups = [v.values for _, v in subdf.groupby('Age')['percent']]
    if len(groups) > 1:  # only if there are multiple response groups
        stat, p = f_oneway(*groups)
        pvals[ct] = p


# Annotate p-values above bars
ax = plt.gca()
ymax = counts['percent'].max()
for i, (ct, p) in enumerate(pvals.items()):
    xloc = i  # bar center for cluster i
    # adjust y slightly above the tallest bar for that cluster
    bar_top = counts[counts['cytetype_annotation_clusters']==ct]['percent'].max()
    ax.text(
        xloc, bar_top + 1,  # 3% higher than max bar for visibility
        f"{p:.2g}",
        ha='center', va='bottom', fontsize=10
    )

plt.title('Celltype Percentage')
plt.xlabel('Celltype')
plt.ylabel('Cell Percentage (%)')
plt.xticks(rotation=90)
plt.legend(title='Age', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

os.chdir(savedir)
os.makedirs("./figures",exist_ok=True)
plt.savefig("./figures/celltype_response_barplot_anova_pvalue.pdf")
plt.show()

os.makedirs("Table",exist_ok=True)
counts.to_csv("./Table/Celltype_Age_barplot_anova_pvalue.txt")

cluster_celltype=pd.DataFrame((adata_full3.obs['clusters'].astype(str) + "_" + adata_full3.obs['cytetype_annotation_clusters'].astype(str)).unique().tolist())

savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/CD4/Day1/"
cluster_celltype.to_csv(os.path.join(savedir,"Table/cluster_celltype.txt"), index = False, header = None)
adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day1_annotation.h5ad"))

#region Day3
import anndata
import scanpy as sc
from cytetype import CyteType
import mygene
import pandas as pd
import os

# ------ Example Scanpy Pipeline ------
#  Skip this step if you already have clusters and marker genes in an AnnData object. 
adata = anndata.read_h5ad("/mnt/data/projects/.immune/Mayo/Activation_TP/processed_data/CD4_D3.h5ad")
# Since the log normalization is already done
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata_full)

adata_full = adata.raw.to_adata()
### Changing the gene name
mg = mygene.MyGeneInfo()
ensg_ids = adata_full.var_names.tolist()
gene_info = mg.querymany(ensg_ids, scopes='ensembl.gene', fields='symbol', species='human')
gene_map = pd.DataFrame(gene_info)[['query', 'symbol']].dropna()

adata_full.var['gene_symbol'] = adata_full.var_names.map(
    dict(zip(gene_map['query'], gene_map['symbol']))
)

adata_full2 = adata_full[:,~adata_full.var['gene_symbol'].isna()].copy() # remove Nan values
adata_full3 = adata_full2[:,~adata_full2.var['gene_symbol'].duplicated()] # remove duplicated
adata_full3.var.index = adata_full3.var['gene_symbol']

adata_full3.obs['clusters'] = adata_full3.obs['seurat_clusters'].astype(str).str.replace("Cluster ","")
sc.tl.rank_genes_groups(adata_full3, groupby="clusters", method="t-test")

adata_full3.var['gene_symbols'] = adata_full3.var.index
annotator = CyteType(adata_full3, group_key="clusters")
adata_full3 = annotator.run(
    study_context="CD4+ T cells extracted from PBMC of healthy individuals stimulated for 3 days with CD3-CD28 beads",
)

savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/CD4/Day3/"
os.chdir(savedir)
sc.pl.umap(adata_full3,color="cytetype_annotation_clusters",size=5, save = "_CD4_Day3_cytetype_annotation.pdf")
sc.pl.umap(adata_full3,color="clusters",size=5, save = "_CD4_Day3_clusters.pdf")
sc.pl.umap(adata_full3,color="seurat_clusters",size=5, save = "_CD4_Day3_seurat_clusters.pdf")
sc.pl.umap(adata_full3,color="samples",size=4, save = "_CD4_Day3_sample.pdf")

# sc.pl.dotplot(adata_full3,markers,'cytetype_annotation_clusters', standard_scale = "var",save = "_cytetype_annotation_2.pdf")
# sc.pl.heatmap(adata_full3,markers,'cytetype_annotation_clusters',standard_scale = "var",save = "_cytetype_annotation_2.pdf")
os.makedirs(os.path.join(savedir,"saveh5ad"),exist_ok = True)
adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day3_annotation.h5ad"))

sc.tl.rank_genes_groups(adata_full3, groupby="cytetype_annotation_clusters", method="t-test", pts=True)

# Extract results
res = adata_full3.uns['rank_genes_groups']
groups = res['names'].dtype.names  # all cell type names

# Define save directory
os.makedirs(os.path.join(savedir,"Table"), exist_ok=True)

# Loop through all cell types
for group in groups:
    genes = res['names'][group]
    # build DataFrame
    de_df = pd.DataFrame({
        'gene': genes,
        'logfoldchange': res['logfoldchanges'][group],
        'pvals': res['pvals'][group],
        'pvals_adj': res['pvals_adj'][group],
    })
    # add expression percentages (if pts=True was used)
    if 'pts' in res:
        pts = pd.DataFrame(res['pts'], index=adata_full3.var_names)
        pts_rest = pd.DataFrame(res['pts_rest'], index=adata_full3.var_names)  
        # subset to genes in this group
        pts_group = pts.loc[genes]
        pts_rest_group = pts_rest.loc[genes]
        de_df['pct_expr_in_group'] = pts_group[group].values
        de_df['pct_expr_in_rest'] = pts_rest_group[group].values
    # save to CSV
    group=group.replace("/","_")
    out_path = os.path.join(savedir,"Table", f"{group}_DE_results.csv")
    de_df.to_csv(out_path, index=False)
    print(f"Saved {group} markers → {out_path}")

#region quantification
import pandas as pd
import os

df = adata_full3.obs[['samples','cytetype_annotation_clusters']].copy()
counts = df.groupby(['samples','cytetype_annotation_clusters']).size().reset_index(name='cell_count')
counts['Age'] = counts['samples'].astype(str).str.replace(r'(\D+)\d.*', r'\1', regex=True)

counts['total_cells'] = counts.groupby('samples')['cell_count'].transform('sum')
counts['percent'] = (counts['cell_count'] / counts['total_cells']) * 100

counts['Age'] = pd.Categorical(counts['Age'],categories=['HY','HO','FO'],ordered=True)

## Just checking
counts.loc[counts['samples'].isin(["HY4D1N4"]),"percent"].sum()

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 8))
sns.barplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    palette='Set2',
    errorbar='sd'
)

sns.stripplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    dodge=True,
    palette='Set2',
    linewidth=0.8,
    edgecolor='black',
    size=5,
    marker='o',
    alpha=0.7,
    legend=False  # keep legend from barplot only
)

from scipy.stats import f_oneway
pvals = {}
for ct, subdf in counts.groupby('cytetype_annotation_clusters'):
    groups = [v.values for _, v in subdf.groupby('Age')['percent']]
    if len(groups) > 1:  # only if there are multiple response groups
        stat, p = f_oneway(*groups)
        pvals[ct] = p


# Annotate p-values above bars
ax = plt.gca()
ymax = counts['percent'].max()
for i, (ct, p) in enumerate(pvals.items()):
    xloc = i  # bar center for cluster i
    # adjust y slightly above the tallest bar for that cluster
    bar_top = counts[counts['cytetype_annotation_clusters']==ct]['percent'].max()
    ax.text(
        xloc, bar_top + 1,  # 3% higher than max bar for visibility
        f"{p:.2g}",
        ha='center', va='bottom', fontsize=10
    )

plt.title('Celltype Percentage')
plt.xlabel('Celltype')
plt.ylabel('Cell Percentage (%)')
plt.xticks(rotation=90)
plt.legend(title='Age', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

os.chdir(savedir)
os.makedirs("./figures",exist_ok=True)
plt.savefig("./figures/celltype_response_barplot_anova_pvalue.pdf")
plt.show()

os.makedirs("Table",exist_ok=True)
counts.to_csv("./Table/Celltype_Age_barplot_anova_pvalue.txt")

cluster_celltype=pd.DataFrame((adata_full3.obs['clusters'].astype(str) + "_" + adata_full3.obs['cytetype_annotation_clusters'].astype(str)).unique().tolist())

cluster_celltype.to_csv(os.path.join(savedir,"Table/cluster_celltype.txt"), index = False, header = None)
adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day3_annotation.h5ad"))

#region Day6
import anndata
import scanpy as sc
from cytetype import CyteType
import mygene
import pandas as pd
import os

# ------ Example Scanpy Pipeline ------
#  Skip this step if you already have clusters and marker genes in an AnnData object. 
adata = anndata.read_h5ad("/mnt/data/projects/.immune/Mayo/Activation_TP/processed_data/CD4_D6.h5ad")
# Since the log normalization is already done
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata_full)

adata_full = adata.raw.to_adata()
### Changing the gene name
mg = mygene.MyGeneInfo()
ensg_ids = adata_full.var_names.tolist()
gene_info = mg.querymany(ensg_ids, scopes='ensembl.gene', fields='symbol', species='human')
gene_map = pd.DataFrame(gene_info)[['query', 'symbol']].dropna()

adata_full.var['gene_symbol'] = adata_full.var_names.map(
    dict(zip(gene_map['query'], gene_map['symbol']))
)

adata_full2 = adata_full[:,~adata_full.var['gene_symbol'].isna()].copy() # remove Nan values
adata_full3 = adata_full2[:,~adata_full2.var['gene_symbol'].duplicated()] # remove duplicated
adata_full3.var.index = adata_full3.var['gene_symbol']

adata_full3.obs['clusters'] = adata_full3.obs['seurat_clusters'].astype(str).str.replace("Cluster ","")
sc.tl.rank_genes_groups(adata_full3, groupby="clusters", method="t-test")

adata_full3.var['gene_symbols'] = adata_full3.var.index
annotator = CyteType(adata_full3, group_key="clusters")
adata_full3 = annotator.run(
    study_context="CD4+ T cells extracted from PBMC of healthy individuals stimulated for 6 days with CD3-CD28 beads",
)

savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/CD4/Day6/"
os.makedirs(savedir,exist_ok=True)
os.chdir(savedir)
sc.pl.umap(adata_full3,color="cytetype_annotation_clusters",size=5, save = "_CD4_Day6_cytetype_annotation.pdf")
sc.pl.umap(adata_full3,color="clusters",size=5, save = "_CD4_Day6_clusters.pdf")
sc.pl.umap(adata_full3,color="seurat_clusters",size=5, save = "_CD4_Day6_seurat_clusters.pdf")
sc.pl.umap(adata_full3,color="samples",size=4, save = "_CD4_Day6_sample.pdf")

# sc.pl.dotplot(adata_full3,markers,'cytetype_annotation_clusters', standard_scale = "var",save = "_cytetype_annotation_2.pdf")
# sc.pl.heatmap(adata_full3,markers,'cytetype_annotation_clusters',standard_scale = "var",save = "_cytetype_annotation_2.pdf")
os.makedirs(os.path.join(savedir,"saveh5ad"),exist_ok = True)
sc.tl.rank_genes_groups(adata_full3, groupby="cytetype_annotation_clusters", method="t-test", pts=True)
adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day6_annotation.h5ad"))

# Extract results
res = adata_full3.uns['rank_genes_groups']
groups = res['names'].dtype.names  # all cell type names

# Define save directory
os.makedirs(os.path.join(savedir,"Table"), exist_ok=True)

# Loop through all cell types
for group in groups:
    genes = res['names'][group]
    # build DataFrame
    de_df = pd.DataFrame({
        'gene': genes,
        'logfoldchange': res['logfoldchanges'][group],
        'pvals': res['pvals'][group],
        'pvals_adj': res['pvals_adj'][group],
    })
    # add expression percentages (if pts=True was used)
    if 'pts' in res:
        pts = pd.DataFrame(res['pts'], index=adata_full3.var_names)
        pts_rest = pd.DataFrame(res['pts_rest'], index=adata_full3.var_names)  
        # subset to genes in this group
        pts_group = pts.loc[genes]
        pts_rest_group = pts_rest.loc[genes]
        de_df['pct_expr_in_group'] = pts_group[group].values
        de_df['pct_expr_in_rest'] = pts_rest_group[group].values
    # save to CSV
    group=group.replace("/","_")
    out_path = os.path.join(savedir,"Table", f"{group}_DE_results.csv")
    de_df.to_csv(out_path, index=False)
    print(f"Saved {group} markers → {out_path}")

#region quantification
import pandas as pd
import os

df = adata_full3.obs[['samples','cytetype_annotation_clusters']].copy()
counts = df.groupby(['samples','cytetype_annotation_clusters']).size().reset_index(name='cell_count')
counts['Age'] = counts['samples'].astype(str).str.replace(r'(\D+)\d.*', r'\1', regex=True)

counts['total_cells'] = counts.groupby('samples')['cell_count'].transform('sum')
counts['percent'] = (counts['cell_count'] / counts['total_cells']) * 100

counts['Age'] = pd.Categorical(counts['Age'],categories=['HY','HO','FO'],ordered=True)

## Just checking
counts.loc[counts['samples'].isin(["HY4D1N4"]),"percent"].sum()

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 8))
sns.barplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    palette='Set2',
    errorbar='sd'
)

sns.stripplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    dodge=True,
    palette='Set2',
    linewidth=0.8,
    edgecolor='black',
    size=5,
    marker='o',
    alpha=0.7,
    legend=False  # keep legend from barplot only
)

from scipy.stats import f_oneway
pvals = {}
for ct, subdf in counts.groupby('cytetype_annotation_clusters'):
    groups = [v.values for _, v in subdf.groupby('Age')['percent']]
    if len(groups) > 1:  # only if there are multiple response groups
        stat, p = f_oneway(*groups)
        pvals[ct] = p


# Annotate p-values above bars
ax = plt.gca()
ymax = counts['percent'].max()
for i, (ct, p) in enumerate(pvals.items()):
    xloc = i  # bar center for cluster i
    # adjust y slightly above the tallest bar for that cluster
    bar_top = counts[counts['cytetype_annotation_clusters']==ct]['percent'].max()
    ax.text(
        xloc, bar_top + 1,  # 3% higher than max bar for visibility
        f"{p:.2g}",
        ha='center', va='bottom', fontsize=10
    )

plt.title('Celltype Percentage')
plt.xlabel('Celltype')
plt.ylabel('Cell Percentage (%)')
plt.xticks(rotation=90)
plt.legend(title='Age', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

os.chdir(savedir)
os.makedirs("./figures",exist_ok=True)
plt.savefig("./figures/celltype_response_barplot_anova_pvalue.pdf")
plt.show()

os.makedirs("Table",exist_ok=True)
counts.to_csv("./Table/Celltype_Age_barplot_anova_pvalue.txt")

cluster_celltype=pd.DataFrame((adata_full3.obs['clusters'].astype(str) + "_" + adata_full3.obs['cytetype_annotation_clusters'].astype(str)).unique().tolist())

cluster_celltype.to_csv(os.path.join(savedir,"Table/cluster_celltype.txt"), index = False, header = None)
# adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day6_annotation.h5ad"))

#region Day0
import anndata
import scanpy as sc
from cytetype import CyteType
import mygene
import pandas as pd
import os

# ------ Example Scanpy Pipeline ------
#  Skip this step if you already have clusters and marker genes in an AnnData object. 
adata = anndata.read_h5ad("/mnt/data/projects/.immune/Mayo/Activation_TP/processed_data/CD4_D0.h5ad")
# Since the log normalization is already done
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata_full)

adata_full = adata.raw.to_adata()
### Changing the gene name
mg = mygene.MyGeneInfo()
ensg_ids = adata_full.var_names.tolist()
gene_info = mg.querymany(ensg_ids, scopes='ensembl.gene', fields='symbol', species='human')
gene_map = pd.DataFrame(gene_info)[['query', 'symbol']].dropna()

adata_full.var['gene_symbol'] = adata_full.var_names.map(
    dict(zip(gene_map['query'], gene_map['symbol']))
)

adata_full2 = adata_full[:,~adata_full.var['gene_symbol'].isna()].copy() # remove Nan values
adata_full3 = adata_full2[:,~adata_full2.var['gene_symbol'].duplicated()] # remove duplicated
adata_full3.var.index = adata_full3.var['gene_symbol']

adata_full3.obs['clusters'] = adata_full3.obs['seurat_clusters'].astype(str).str.replace("Cluster ","")
sc.tl.rank_genes_groups(adata_full3, groupby="clusters", method="t-test")

adata_full3.var['gene_symbols'] = adata_full3.var.index
annotator = CyteType(adata_full3, group_key="clusters")
adata_full3 = annotator.run(
    study_context="CD4+ T cells extracted from PBMC of healthy individuals at the resting phase",
)

savedir = "/mnt/data/projects/.immune/Mayo/Activation_TP/analysis/CD4/Day0/"
os.makedirs(savedir,exist_ok=True)
os.chdir(savedir)
sc.pl.umap(adata_full3,color="cytetype_annotation_clusters",size=5, save = "_CD4_Day0_cytetype_annotation.pdf")
sc.pl.umap(adata_full3,color="clusters",size=5, save = "_CD4_Day0_clusters.pdf")
sc.pl.umap(adata_full3,color="seurat_clusters",size=5, save = "_CD4_Day0_seurat_clusters.pdf")
sc.pl.umap(adata_full3,color="samples",size=4, save = "_CD4_Day0_sample.pdf")

# sc.pl.dotplot(adata_full3,markers,'cytetype_annotation_clusters', standard_scale = "var",save = "_cytetype_annotation_2.pdf")
# sc.pl.heatmap(adata_full3,markers,'cytetype_annotation_clusters',standard_scale = "var",save = "_cytetype_annotation_2.pdf")
os.makedirs(os.path.join(savedir,"saveh5ad"),exist_ok = True)
adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day0_annotation.h5ad"))

sc.tl.rank_genes_groups(adata_full3, groupby="cytetype_annotation_clusters", method="t-test", pts=True)
# Extract results
res = adata_full3.uns['rank_genes_groups']
groups = res['names'].dtype.names  # all cell type names

# Define save directory
os.makedirs(os.path.join(savedir,"Table"), exist_ok=True)

# Loop through all cell types
for group in groups:
    genes = res['names'][group]
    # build DataFrame
    de_df = pd.DataFrame({
        'gene': genes,
        'logfoldchange': res['logfoldchanges'][group],
        'pvals': res['pvals'][group],
        'pvals_adj': res['pvals_adj'][group],
    })
    # add expression percentages (if pts=True was used)
    if 'pts' in res:
        pts = pd.DataFrame(res['pts'], index=adata_full3.var_names)
        pts_rest = pd.DataFrame(res['pts_rest'], index=adata_full3.var_names)  
        # subset to genes in this group
        pts_group = pts.loc[genes]
        pts_rest_group = pts_rest.loc[genes]
        de_df['pct_expr_in_group'] = pts_group[group].values
        de_df['pct_expr_in_rest'] = pts_rest_group[group].values
    # save to CSV
    group=group.replace("/","_")
    out_path = os.path.join(savedir,"Table", f"{group}_DE_results.csv")
    de_df.to_csv(out_path, index=False)
    print(f"Saved {group} markers → {out_path}")

#region quantification
import pandas as pd
import os

df = adata_full3.obs[['samples','cytetype_annotation_clusters']].copy()
counts = df.groupby(['samples','cytetype_annotation_clusters']).size().reset_index(name='cell_count')
counts['Age'] = counts['samples'].astype(str).str.replace(r'(\D+)\d.*', r'\1', regex=True)

counts['total_cells'] = counts.groupby('samples')['cell_count'].transform('sum')
counts['percent'] = (counts['cell_count'] / counts['total_cells']) * 100

counts['Age'] = pd.Categorical(counts['Age'],categories=['HY','HO','FO'],ordered=True)

## Just checking
counts.loc[counts['samples'].isin(["HY4D1N4"]),"percent"].sum()

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 8))
sns.barplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    palette='Set2',
    errorbar='sd'
)

sns.stripplot(
    data=counts,
    x='cytetype_annotation_clusters',
    y='percent',
    hue='Age',
    dodge=True,
    palette='Set2',
    linewidth=0.8,
    edgecolor='black',
    size=5,
    marker='o',
    alpha=0.7,
    legend=False  # keep legend from barplot only
)

from scipy.stats import f_oneway
pvals = {}
for ct, subdf in counts.groupby('cytetype_annotation_clusters'):
    groups = [v.values for _, v in subdf.groupby('Age')['percent']]
    if len(groups) > 1:  # only if there are multiple response groups
        stat, p = f_oneway(*groups)
        pvals[ct] = p


# Annotate p-values above bars
ax = plt.gca()
ymax = counts['percent'].max()
for i, (ct, p) in enumerate(pvals.items()):
    xloc = i  # bar center for cluster i
    # adjust y slightly above the tallest bar for that cluster
    bar_top = counts[counts['cytetype_annotation_clusters']==ct]['percent'].max()
    ax.text(
        xloc, bar_top + 1,  # 3% higher than max bar for visibility
        f"{p:.2g}",
        ha='center', va='bottom', fontsize=10
    )

plt.title('Celltype Percentage')
plt.xlabel('Celltype')
plt.ylabel('Cell Percentage (%)')
plt.xticks(rotation=90)
plt.legend(title='Age', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

os.chdir(savedir)
os.makedirs("./figures",exist_ok=True)
plt.savefig("./figures/celltype_response_barplot_anova_pvalue.pdf")
plt.show()

os.makedirs("Table",exist_ok=True)
counts.to_csv("./Table/Celltype_Age_barplot_anova_pvalue.txt")

cluster_celltype=pd.DataFrame((adata_full3.obs['clusters'].astype(str) + "_" + adata_full3.obs['cytetype_annotation_clusters'].astype(str)).unique().tolist())

cluster_celltype.to_csv(os.path.join(savedir,"Table/cluster_celltype.txt"), index = False, header = None)
# adata_full3.write_h5ad(os.path.join(savedir,"saveh5ad/CD4_Day0_annotation.h5ad"))
