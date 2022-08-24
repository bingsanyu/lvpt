from multiprocessing import freeze_support
import scvelo as scv
import scanpy as sc
import csv
import lvpt
from velocity_graph import velocity_graph

def main():
    scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
    scv.settings.presenter_view = True  # set max width size for presenter view
    scv.set_figure_params('scvelo')  # for beautified visualization

    # Read simulated data
    adata = scv.read("data/linear.h5ad", cache=True)
    # adata = scv.read("data/bifurcating.h5ad", cache=True)
    # adata = scv.read("data/trifurcating.h5ad", cache=True)
    adata.layers['spliced']=adata.layers['counts_spliced']
    adata.layers['unspliced']=adata.layers['counts_unspliced']

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.leiden(adata, resolution=0.54)
    
    # Velocity
    scv.tl.velocity(adata)
    # Pseudotime
    velocity_graph(adata, steady = 0.05)
    # linear
    lvpt.velocity_pseudotime(adata, root_key='cell9')
    # bifurcating
    # lvpt.velocity_pseudotime(adata, root_key='cell68')
    # trifurcating
    # lvpt.velocity_pseudotime(adata, root_key='cell858')

    # Save Pseudotime
    with open ('simulated/lvpt-linear.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerow(adata.obs['velocity_pseudotime'])


if __name__ == '__main__':
    freeze_support()
    main()