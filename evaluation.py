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
    nums = [1,11,21,61,81]
    method = "lvpt"
    for i in nums:
        if i%10==1:
            prefix = "linear"
        elif i%10==2:
            prefix = "bifurcating"
        elif i%10==3 :
            prefix = "trifurcating"
        filename = prefix+str(i)
            
        adata = scv.read("simulated_data/"+filename+".h5ad", cache=True)
        adata.layers['spliced']=adata.layers['counts_spliced']
        adata.layers['unspliced']=adata.layers['counts_unspliced']
        
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        sc.tl.leiden(adata, resolution=0.54)
        
        # Calculate velocity
        scv.tl.velocity(adata)
        # Calculate pseudotime
        if(method=="lvpt"):
            velocity_graph(adata, steady = 0.05)
        elif(method=="vpt"):
            scv.tl.velocity_graph(adata)
        scv.tl.terminal_states(adata)
        root = adata.obs['root_cells'].idxmax()
        with open("simulated_result/root_"+str(i)+".txt", 'w') as f:
            f.write(root)

        if(method=="lvpt"):
            lvpt.velocity_pseudotime(adata, root_key=root)
        elif(method=="vpt"):
            scv.tl.velocity_pseudotime(adata, root_key=root)

        filename = filename+method
        # Save pseudotime
        with open ("simulated_result/"+filename+".csv",'w',newline = '') as csvfile:
            my_writer = csv.writer(csvfile, delimiter = ',')
            my_writer.writerow(adata.obs['velocity_pseudotime'])

if __name__ == '__main__':
    freeze_support()
    main()