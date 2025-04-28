# generating shuffled cytokine networks for validation
# the resulting shuffled objects are re-ran with the pipeline, just replacing the input rds files

library(tidyverse)
library(Seurat)
library(nichenetr)

source('iteratecells.R')

# Read nichenet model
ligand_target_matrix <- readRDS("../rds/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("../rds/weighted_networks_nsga2r_final.rds")
ligand_tf_matrix <- readRDS('../rds/ligand_tf_matrix_nsga2r_final.rds')
networks <- readRDS("../rds/OP_nichenet_networks.RDS")
lr_network <- networks$lr_network %>% dplyr::filter(from %in% colnames(ligand_target_matrix))
lr_network <- dplyr::rename(lr_network, ligand = from, receptor = to)
sig_network <- networks$signaling_network
gr_network <- networks$gr_network
geneset_oi <- read_tsv("../rds/cytokine_list.tsv") %>%
  pull(Name) %>%
  unique() %>%
  .[. %in% rownames(ligand_target_matrix)] %>%
  .[. %in% colnames(ligand_target_matrix)]

# Function to rename the integrated single-cell data, function reused from https://github.com/vertesy/Seurat.utils
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { 
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

# setting seed for reproducibility
set.seed(42)

# UC
# inflamed_naive

uc <- readRDS('../rds/uc_inflamed_naive.rds')
uc<- NormalizeData(uc)
Idents(uc) <-'minor_cluster'
# original gene names from expression matrix
og_gene_names <- rownames(uc@assays$RNA$counts)
# shuffle character vector using sample
shuffled_gene_names<-sample(og_gene_names)
# new seurat object with shuffled gene names
uc2<-RenameGenesSeurat(obj = uc, newnames = shuffled_gene_names)
rm(uc)
celltypes = uc2 %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc2)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_inflamed_naive_scibd_shuffled.tsv')

# inflamed treated
uc <- readRDS('../rds/uc_inflamed_no_naive.rds')
uc<- NormalizeData(uc)
Idents(uc) <-'minor_cluster'
# original gene names from expression matrix
og_gene_names <- rownames(uc@assays$RNA$counts)
# shuffle character vector using sample
shuffled_gene_names<-sample(og_gene_names)
# new seurat object with shuffled gene names
uc2<-RenameGenesSeurat(obj = uc, newnames = shuffled_gene_names)
rm(uc)
uc2<- NormalizeData(uc2)
celltypes = uc2 %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc2)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_inflamed_treated_scibd_shuffled.tsv')


# noninflamed treated
uc <- readRDS('../rds/uc_noninflamed_no_naive.rds')
uc<- NormalizeData(uc)
Idents(uc) <-'minor_cluster'
# original gene names from expression matrix
og_gene_names <- rownames(uc@assays$RNA$counts)
# shuffle character vector using sample
shuffled_gene_names<-sample(og_gene_names)
# new seurat object with shuffled gene names
uc2<-RenameGenesSeurat(obj = uc, newnames = shuffled_gene_names)
rm(uc)
uc2<- NormalizeData(uc2)
celltypes = uc2 %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc2)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'uc_noninflamed_treated_scibd_shuffled.tsv')

# Healthy
uc <- readRDS('../rds/scibd_healthy.rds')
uc<- NormalizeData(uc)
Idents(uc) <-'minor_cluster'
# original gene names from expression matrix
og_gene_names <- rownames(uc@assays$RNA$counts)
# shuffle character vector using sample
shuffled_gene_names<-sample(og_gene_names)
# new seurat object with shuffled gene names
uc2<-RenameGenesSeurat(obj = uc, newnames = shuffled_gene_names)
rm(uc)
uc2<- NormalizeData(uc2)
celltypes = uc2 %>% pull(minor_cluster) |> as.character() |> unique()

results <- list()
for (i in celltypes){
  for (j in celltypes){
    hits<-iterateCells(i,j, uc2)
    #print(hits)
    results[[paste0(i,":",j)]]<-hits
  }
}
ucres<-keep(results, is.data.frame) |> bind_rows()
write_tsv(ucres,'healthy_scibd_shuffled.tsv')
