# Generating all CCI models using NicheNet
# Author: Marton Olbei

setwd('~/Documents/Projects/CKL2/nichenet/scIBD/scripts/')
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tidyseurat)
library(nichenetr)
library(ggraph)
library(tidygraph)
library(igraph)
source('iteratecells.R')

# nichenet
ligand_target_matrix <- readRDS("../rds/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("../rds/weighted_networks_nsga2r_final.rds")
ligand_tf_matrix <- readRDS('../rds/ligand_tf_matrix_nsga2r_final.rds')
networks <- readRDS("../rds/OP_nichenet_networks.RDS")
lr_network <- networks$lr_network %>% dplyr::filter(from %in% colnames(ligand_target_matrix))
lr_network <- dplyr::rename(lr_network, ligand = from, receptor = to)
sig_network <- networks$signaling_network
gr_network <- networks$gr_network
geneset_oi <- read_tsv("../cytokine_list.tsv") %>%
  pull(Name) %>%
  unique() %>%
  .[. %in% rownames(ligand_target_matrix)] %>%
  .[. %in% colnames(ligand_target_matrix)]


# UC
# inflamed_naive
uc <- readRDS('../rds/uc_inflamed_naive.rds')
Idents(uc) <-'minor_cluster'

uc<- NormalizeData(uc)

celltypes = uc %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_inflamed_scibd_p01_weights_naive_new_model_input.tsv')

# inflamed_no_naive
uc <- readRDS('../rds/uc_inflamed_no_naive.rds')
Idents(uc) <-'minor_cluster'

uc<- NormalizeData(uc)

celltypes = uc %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_inflamed_scibd_p01_weights_no_naive_new_model_input.tsv')

# noninflamed_naive
uc <- readRDS('../rds/uc_noninflamed_naive.rds')
Idents(uc) <-'minor_cluster'

uc<- NormalizeData(uc)

celltypes = uc %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_noninflamed_scibd_p01_weights_naive_new_model_input.tsv')

# noninflamed_no_naive
uc <- readRDS('../rds/uc_noninflamed_no_naive.rds')
Idents(uc) <-'minor_cluster'

uc<- NormalizeData(uc)

celltypes = uc %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_noninflamed_scibd_p01_weights_no_naive_new_model_input.tsv')

# Healthy
uc <- readRDS('../rds/scibd_healthy.rds')
Idents(uc) <-'minor_cluster'

uc<- NormalizeData(uc)

celltypes = uc %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'healthy_scibd_p01_weights_new_model_input.tsv')

