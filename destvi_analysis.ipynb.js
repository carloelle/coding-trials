import destvi_utils
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy
from scvi.model import CondSCVI, DestVI
from skmisc.loess import loess
import torch
    

sc_adata=sc.read_h5ad('scRNA_mouse_PDAC_day30.h5ad')
st_adata=sc.read_visium('pathtofile/outs',source_image_path='pathtofiletv/outs/spatial')
st_filtered=pd.read_csv('SelectedSpots.csv')


# subset st data
st_adata.var_names_make_unique()
st_adata=st_adata[st_filtered['x'],]
st_adata

# NB: sc_adata contains raw counts
sc.pp.filter_genes(sc_adata, min_counts=10)
G = 2000
sc_adata.layers["counts"] = sc_adata.X.copy()
sc.pp.highly_variable_genes(sc_adata, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3")
sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

# Spatial data
st_adata.layers["counts"] = st_adata.X.copy()

sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

loc=st_adata.obsm["spatial"]
st_adata.obsm["spatial"]=loc.astype('float')

# filter genes to be the same on the spatial and sc data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()

# Fit the scLMV
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="Annotation")
sc_model = CondSCVI(sc_adata, weight_obs=False)
sc_model.view_anndata_setup()
sc_model.train()

sc_model.history["elbo_train"].iloc[5:].plot()
plt.show()

# Deconvolution
DestVI.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.view_anndata_setup()
st_model.train(max_epochs=2500)
st_model.history["elbo_train"].iloc[10:].plot()
plt.show()

# Get proportions
st_adata.obsm["proportions"] = st_model.get_proportions()
st_adata.obsm["proportions"].to_csv('CellProp_DestVI.csv')

ct_thresholds = destvi_utils.automatic_proportion_threshold(st_adata,  kind_threshold="primary")
ct_thresholds['MonoMacro'] = 0.1

for ct, g in st_model.get_gamma().items():
    st_adata.obsm[f"{ct}_gamma"] = g

# LOAD FUNCTIONS FROM destvi_utils
import anndata as ad
import hotspot
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splev, splrep
from scipy.spatial.distance import pdist, squareform
from sklearn.mixture import GaussianMixture


def _prettify_axis(ax, spatial=False):
    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    if spatial:
        plt.xticks([])
        plt.yticks([])
        plt.xlabel("Spatial1")
        plt.ylabel("Spatial2")


def _form_stacked_quantiles(data, N=100):
    quantiles = np.quantile(data, np.linspace(0, 1, N, endpoint=False))
    return quantiles, np.vstack([_flatten(data, q) for q in quantiles])


def _flatten(x, threshold):
    return (x > threshold) * x


def _smooth_get_critical_points(x, noisy_data, k=5, s=0.1):
    f = splrep(x, noisy_data, k=5, s=1)
    smoothed = splev(x, f)
    derivative = splev(x, f, der=1)
    sign_2nd = splev(x, f, der=2) > 0
    curvature = splev(x, f, der=3)
    return noisy_data, smoothed, derivative, sign_2nd, curvature


def _get_autocorrelations(st_adata, stacked_quantiles, quantiles):
    # create Anndata and run hotspot
    adata = ad.AnnData(stacked_quantiles.T)
    adata.obs_names = st_adata.obs.index
    adata.var_names = [str(i) for i in quantiles]
    adata.obsm["spatial"] = st_adata.obsm["spatial"]
    hs = hotspot.Hotspot(adata, model="none", latent_obsm_key="spatial")
    hs.create_knn_graph(
        weighted_graph=True,
        n_neighbors=10,
    )
    hs_results = hs.compute_autocorrelations(jobs=1)
    index = np.array([float(i) for i in hs_results.index.values])
    return index, hs_results["Z"].values


def _get_laplacian(s, pi):
    N = s.shape[0]
    dist_table = pdist(s)
    bandwidth = np.median(dist_table)
    sigma = 0.5 * bandwidth**2

    l2_square = squareform(dist_table) ** 2
    D = np.exp(-l2_square / sigma) * np.dot(pi, pi.T)
    L = -D
    sum_D = np.sum(D, axis=1)
    for i in range(N):
        L[i, i] = sum_D[i]
    return L


def _get_spatial_components(locations, proportions, data):
    # find top two spatial principal vectors
    # form laplacian
    L = _get_laplacian(locations, proportions)
    # center data
    transla_ = data.copy()
    transla_ -= np.mean(transla_, axis=0)
    # get eigenvectors
    A = np.dot(transla_.T, np.dot(L, transla_))
    w, v = np.linalg.eig(A)
    # don't forget to sort them...
    idx = np.argsort(w)[::-1]
    vec = v[:, idx][:, :]
    return vec


def _vcorrcoef(X, y):
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum((X - Xm) * (y - ym), axis=1)
    r_den = np.sqrt(np.sum((X - Xm) ** 2, axis=1) * np.sum((y - ym) ** 2))
    r = np.divide(
        r_num,
        r_den,
        out=np.zeros_like(
            r_num,
        ),
        where=r_den != 0,
    )
    return r


def _get_delta(lfc):
    return np.max(
        np.abs(GaussianMixture(n_components=3).fit(np.array(lfc).reshape(-1, 1)).means_)
    )

# Get 5 Spatial PCs
gamma = st_model.get_gamma(return_numpy=True)
filter_ = st_adata.obsm["proportions"]['MonoMacro'].values > ct_thresholds['MonoMacro']
locations = st_adata.obsm["spatial"][filter_]
proportions = st_adata.obsm["proportions"]['MonoMacro'].values[filter_]
ct_index = np.where('MonoMacro' == st_model.cell_type_mapping)[0][0]
data = gamma[:, :, ct_index][filter_]

vec=get_spatial_components(locations, proportions, data)[:,:]
projection = np.dot(data - np.mean(data, 0), vec)

SpatialPCs=pd.DataFrame(projection)
SpatialPCs.index=st_adata.obs_names[filter_]
SpatialPCs.to_csv('SpatialPCs_MonoMacro.csv')

# Get genes whose expression correlates with Spatial PCs
sc_adata_slice = sc_adata[sc_adata.obs["Annotation"] == 'MonoMacro']
is_sparse = scipy.sparse.issparse(sc_adata_slice.X)
normalized_counts = sc_adata_slice.X.A if is_sparse else sc_adata_slice.X

indices_ct = np.where(sc_adata.obs["Annotation"] == 'MonoMacro')[0]
sc_latent = sc_model.get_latent_representation(indices=indices_ct)
sc_projection = np.dot(sc_latent - np.mean(sc_latent,0), vec)

r = _vcorrcoef(normalized_counts.T, sc_projection[:, 0])
ranking = np.argsort(r)
PC1Pos=pd.DataFrame(r[ranking][::-1][:50])
PC1Pos.index=list(st_adata.var.index[ranking[::-1][:50]])

PC1Neg=pd.DataFrame(r[ranking][:50])
PC1Neg.index=list(st_adata.var.index[ranking[:50]])

## Generate Expression Matrix for CellType

# impute 
imp_ge = st_model.get_scale_for_ct("MonoMacro", indices=np.where(filter_)[0]).values

# get statistics
avg_library_size = np.mean(np.sum(st_adata.layers["counts"], axis=1).A.flatten())
exp_px_o = st_model.module.px_o.detach().exp().cpu().numpy()
mean = avg_library_size * imp_ge

# create distribution
concentration = torch.tensor(avg_library_size * imp_ge / exp_px_o)
rate = torch.tensor(1. / exp_px_o)

# generate
for j in [1,2,3,4,5,6]:
    N = 1
    simulated = torch.distributions.Gamma(concentration=concentration, rate = rate).sample((N,)).cpu().numpy()
    simulated = np.log(simulated + 1)
    simulated = simulated.reshape((-1, simulated.shape[-1]))
    simulated=pd.DataFrame(simulated, index=st_adata.obs['_indices'][np.where(filter_)[0]].index, columns=st_adata.var['gene_ids'].index)
    simulated.to_csv(f"Simulation_{j}_MonoMacro.csv")