setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dyno)
library(dyneval)
library(dynplot)
library(tidyverse)
library(dyngen)
library(anndata)
library(Matrix)

#Create simulation models

set.seed(1)
backbone <- backbone_linear()

set.seed(2)
backbone <- backbone_bifurcating()

set.seed(6)
backbone <- backbone_binary_tree(
  num_modifications = 3
)

set.seed(4)
backbone <- backbone_cycle()

set.seed(5)
backbone <- backbone_disconnected()

set.seed(9)
backbone <- backbone_trifurcating()

config <- initialise_model(
  backbone = backbone,
  num_cells = 1000,
  num_tfs = 100,
  num_targets = 0,
  num_hks = 0,
  simulation_params = simulation_default(census_interval = 10, ssa_algorithm = ssa_etl(tau = 300 / 3600)),
  verbose = FALSE
)

model <- generate_tf_network(config)

model <- generate_feature_network(model)

model <- generate_kinetics(model)

model <- generate_gold_standard(model)

model <- generate_cells(model)

model <- generate_experiment(model)

dataset <- as_dyno(model)

# Save models
linear = model
save(linear, file = "simulated/linear.RData")

bifurcating = model
save(bifurcating, file = "simulated/bifurcating.RData")

binaryTree = model
save(binaryTree, file = "simulated/binaryTree.RData")

cycle = model
save(cycle, file="simulated/cycle.RData")

disconnected = model
save(disconnected, file="simulated/disconnected.RData")

trifurcating = model
save(trifurcating, file="simulated/trifucating.RData")

# Store AnnData
ad <- as_anndata(bifurcating)
ad$write_h5ad("bifurcating.h5ad")

# Evaluation

load("simulated/linear.RData")
model = linear
dataset = as_dyno(model)
dataset[["prior_information"]][["start_id"]]="cell9"

load("simulated/bifurcating.RData")
model = bifurcating
dataset = as_dyno(model)
dataset[["prior_information"]][["start_id"]]="cell68"

load("simulated/trifucating.RData")
model = trifurcating
dataset = as_dyno(model)
dataset[["prior_information"]][["start_id"]]="cell858"

dataset$pseudotime <- dynwrap::calculate_pseudotime(dataset)

models <- infer_trajectories(dataset, list(ti_monocle_ddrtree(),ti_slingshot(), ti_dpt(), ti_tscan()))
dataset <- add_cell_waypoints(dataset)
models$model <- map(models$model, add_cell_waypoints)
metric_ids <- dyneval::metrics %>% filter(category != "average") %>% pull(metric_id)
metric_ids = metric_ids[-9]
metrics <- map_dfr(models$model, dyneval::calculate_metrics, dataset = dataset, metrics = metric_ids)


# Read lvpt result
time = read.csv("simulated/lvpt-linear.csv",header = F)
time = read.csv("simulated/lvpt-bifurcating.csv",header = F)
time = read.csv("simulated/lvpt-trifurcating.csv",header = F)

time = read.csv("simulated/vpt-linear.csv",header = F)
time = read.csv("simulated/vpt-bifurcating.csv",header = F)
time = read.csv("simulated/vpt-trifurcating.csv",header = F)

cellname = c()
for (i in 1:1000) {
  cellname = c(cellname,paste('cell',i,sep = ''))
}
pseudotime = c()
for(i in time){
  pseudotime = c(pseudotime, i)
}
names(pseudotime) = cellname

cor(pseudotime, dataset$pseudotime, method="spearman")

lvpt = dataset
lvpt$pseudotime = pseudotime

plot_dimred(lvpt,dimred = dyndimred::dimred_pca(dataset$expression), color_cells = "pseudotime",plot_trajector=FALSE)
