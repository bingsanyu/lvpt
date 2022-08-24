from multiprocessing import freeze_support
import scvelo as scv
import scanpy as sc
import lvpt
from velocity_graph import velocity_graph

def main():
    # Preprocessing
    scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
    scv.settings.presenter_view = True  # set max width size for presenter view
    scv.set_figure_params('scvelo')  # for beautified visualization
    adata = scv.datasets.dentategyrus_lamanno()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=100)

    # Clustering
    sc.tl.leiden(adata, resolution=0.85)# 100, 0.85
    colors = ['#a7ff83','#769fcd','#15b7b9','#f38181','#ff6f3c','#95e1d3','#fcbad3','#fce38a','#ffaa64',
                '#17b978','#aa96da',  '#f2c6b4', '#f9a1bc']
    scv.pl.scatter(adata, color='leiden', palette=colors)

    # Velocity
    scv.tl.velocity(adata)
    velocity_graph(adata, steady=0.05, approx=True)
    
    lvpt.velocity_pseudotime(adata, root_key='10X83_2:AAAGTAGCATAAAGGTx')
    scv.pl.scatter(adata, basis='tsne', color='velocity_pseudotime', cmap='rainbow')

    # Differential genes
    scv.tl.rank_velocity_genes(adata, groupby='leiden', min_corr=.3)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.head()
    scv.pl.velocity(adata, ['Gng12'])
    scv.pl.velocity(adata, ['Sema5a'])
    scv.pl.velocity(adata, ['Sema3c'])

    # Trajectory
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    sc.tl.paga(adata, groups='leiden')

    clustertime=adata.obs.groupby("leiden").agg("mean")['velocity_pseudotime']
    for i in adata.obs['leiden'].cat.categories:
        for j in adata.obs['leiden'].cat.categories:
            if(clustertime[i]>clustertime[j]):
                adata.uns['paga']['connectivities'][int(j),int(i)]=0

    scv.pl.paga(adata, basis='tsne', threshold=0.2, transitions='connectivities', solid_edges=None, dashed_edges=None, palette=colors)


if __name__ == '__main__':
    freeze_support()
    main()