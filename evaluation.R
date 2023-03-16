setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dyno)
library(dyneval)
library(dynplot)
library(tidyverse)
library(dyngen)
library(anndata)
library(Matrix)

# Create simulation models
# Set seeds [1, 11, 21, 31, 41, 51] for linear
# Set seeds [2, 12, 22, 32, 42, 52] for bifurcating
# Set seeds [3, 13, 23, 33, 43, 53] for trifurcating
genData <- function(seed){
  set.seed(seed)
  if(seed%%10==1){
    prefix = "linear"
    backbone <- backbone_linear()
  }else if(seed%%10==2){
    prefix = "bifurcating"
    backbone <- backbone_bifurcating()
  }else if(seed%%10==3){
    prefix = "trifurcating"
    backbone <- backbone_trifurcating()
  }
  
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
  
  # Save AnnData
  ad <- as_anndata(model)
  ad$write_h5ad(paste("simulated_data/",prefix, seed, ".h5ad", sep = ''))
  
  # Save Rdata
  save(model, file = paste("simulated_data/",prefix,"_", seed,".RData", sep = ''))
}

seeds = list(1,11,21,61,81)
for(seed in seeds){
  genData(seed)
}

# Evaluation
res = c()
method = "lvpt"
for(seed in seeds){
  if(seed%%10==1){
    prefix = "linear"
  }else if(seed%%10==2){
    prefix = "bifurcating"
  }else if(seed%%10==3){
    prefix = "trifurcating"
  }
  load(paste("simulated_data/",prefix,"_",seed,".RData", sep=""))
  dataset = as_dyno(model)
  root = scan(paste("simulated_result/root_",seed,".txt", sep = ''), what = "", sep = "", na.strings = "NA")
  dataset[["prior_information"]][["start_id"]]=root
  dataset = add_root(dataset, root_cell_id = root)
  dataset$pseudotime <- dynwrap::calculate_pseudotime(dataset)
  # plot_dimred(dataset,dimred = dyndimred::dimred_pca(dataset$expression), color_cells = "pseudotime", plot_trajectory = FALSE)
  
  time = read.csv(paste("simulated_result/",prefix,seed,method,".csv", sep = ''), header = F)
  
  cellname = c()
  for (i in 1:length(time)) {
    cellname = c(cellname,paste('cell',i,sep = ''))
  }
  pseudotime = c()
  for(i in time){
    pseudotime = c(pseudotime, i)
  }
  names(pseudotime) = cellname
  
  p = cor(pseudotime, dataset$pseudotime, method="spearman")
  
  res = c(res,p)
}

monocle = c()
slingshot = c()
tscan = c()

seeds = list(1,11,21,61,81)
for(seed in seeds){
  if(seed%%10==1){
    prefix = "linear"
  }else if(seed%%10==2){
    prefix = "bifurcating"
  }else if(seed%%10==3){
    prefix = "trifurcating"
  }
  load(paste("simulated_data/",prefix,"_",seed,".RData", sep=""))
  dataset = as_dyno(model)
  root = scan(paste("simulated_result/root_",seed,".txt", sep = ''), what = "", sep = "", na.strings = "NA")
  dataset[["prior_information"]][["start_id"]]=root
  dataset = add_root(dataset, root_cell_id = root)
  models <- infer_trajectories(dataset, list(ti_monocle_ddrtree(),ti_slingshot(), ti_tscan()))
  dataset <- add_cell_waypoints(dataset)
  models$model <- map(models$model, add_cell_waypoints)
  metric_ids <- dyneval::metrics %>% filter(category != "average") %>% pull(metric_id)
  metric_ids = metric_ids[-9]
  metrics <- map_dfr(models$model, dyneval::calculate_metrics, dataset = dataset, metrics = metric_ids)
  monocle = c(monocle, metrics$correlation[1])
  slingshot = c(slingshot, metrics$correlation[2])
  tscan = c(tscan, metrics$correlation[3])
}

