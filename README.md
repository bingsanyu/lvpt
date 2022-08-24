# LVPT
Lazy velocity pseudotime method for single cell RNA-seq data.

# Prepare
To run lvpt, please install following softwares:

* scanpy
* scvelo
* dynverse
* docker (To make sure the dynmethod can run)

# Start
Use following code to calculate laze velocity pseudotime.

```velocity_graph(adata, steady = 0.05)\nlvpt.velocity_pseudotime(adata, root_key='cell9')```

