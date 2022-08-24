from multiprocessing import freeze_support
import scvelo as scv
import scanpy as sc
import lvpt
from velocity_graph import velocity_graph

def main():
    #Preprocessing
    scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
    scv.settings.presenter_view = True  # set max width size for presenter view
    scv.set_figure_params('scvelo')  # for beautified visualization
    adata = scv.datasets.pancreas()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    #Clustering
    sc.tl.leiden(adata, resolution=0.54)
    colors = ['#a7ff83','#769fcd','#ffaa64','#aa96da','#f38181','#95e1d3','#fcbad3','#fce38a',
                 '#17b978', '#ff6f3c','#f2c6b4','#15b7b9', '#f9a1bc']
    scv.pl.scatter(adata, color='leiden', palette=colors)
    
    # Velocity
    scv.tl.velocity(adata)
    
    # Pseudotime
    velocity_graph(adata, steady = 0.05)
    lvpt.velocity_pseudotime(adata, root_key='AACACGTTCCGCGCAA', end_key='AAAGATGCAATGTTGC')
    scv.pl.scatter(adata, color='velocity_pseudotime', cmap='rainbow')

    # Differential genes
    scv.tl.rank_velocity_genes(adata, groupby='leiden', min_corr=.3)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.head()
    scv.pl.velocity(adata, ['Abcc8'])
    scv.pl.velocity(adata, ['Gnas'])
    
    # Heatmap
    scv.pl.heatmap(adata, ['Notch2','Fgf12','Chd7','Adgrb3','Sdk1','Heg1','Pax6','Zcchc16','Spock3'],sortby='velocity_pseudotime', sort=False)

    # Trajectory
    scv.tl.velocity_graph(adata)
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    scv.tl.paga(adata, groups='leiden')
    scv.pl.paga(adata, basis='umap', solid_edges='connectivities', dashed_edges=None)

    adata

if __name__ == '__main__':
    freeze_support()
    main()