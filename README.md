# LVPT
Lazy velocity pseudotime method for single cell RNA-seq data.

# Prepare
To run lvpt, please install following softwares:

1. scanpy
2. scvelo
3. dynverse
4. docker (To make sure the dynmethod can run)

# Start
Use following code to calculate laze velocity pseudotime.
`velocity_graph(adata, steady = 0.05)`
`lvpt.velocity_pseudotime(adata, root_key='cell9')`

