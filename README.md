# LVPT
Lazy velocity pseudotime method for single cell RNA-seq data.

# Prepare
To run lvpt, please install following softwares:

* scanpy
* scvelo
* dynverse
* docker (To make sure the dynmethod can run)

# Start
Import Lvpt files.

```import lvpt```
```from velocity_graph import velocity_graph```

Use following code to calculate lazy velocity pseudotime.

```velocity_graph(adata, steady = 0.05)```
```lvpt.velocity_pseudotime(adata, root_key='cell9')```

# Examples
File name | Description
--------- | -----------
main.py   | The example experiment of pancreas dataset.
dentategyus.py | The example experiment of dentategyus dataset.
evaluation.py & evaluation.R | The examples experiment of simulated datasets.
