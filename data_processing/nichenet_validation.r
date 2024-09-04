library(nichenetr)
library(tidyverse)
library(Seurat)

#nichenet
ligand_target_matrix <- readRDS("../rds/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("../rds/weighted_networks_nsga2r_final.rds")
ligand_tf_matrix <- readRDS("../rds/ligand_tf_matrix_nsga2r_final.rds")
networks <- readRDS('../rds/OP_nichenet_networks.RDS')
lr_network <- networks$lr_network %>% dplyr::filter(from %in% colnames(ligand_target_matrix))
lr_network<- dplyr::rename(lr_network, ligand = from, receptor = to)
sig_network <- networks$signaling_network
gr_network <- networks$gr_network
## geneset of interest == cytokines
geneset_oi <- read_tsv("../rds/cytokine_list.tsv") %>%
  pull(Name) %>%
  unique()  %>%
  .[. %in% rownames(ligand_target_matrix)]  %>% .[. %in% colnames(ligand_target_matrix)]

# processing function
rf_validation_CKL <- function(rds_file, tsv_file, output_file){
  ucint <- readRDS(paste0('../rds/',rds_file))
  Idents(ucint) <-'minor_cluster'
  ligand_df <- read_tsv(paste0(tsv_file)) %>%  dplyr::select(target_cell,ligand) %>%  distinct()
  tryCatch(
    expr = {

    }, error = function(e) {
      # Code to run when there's an error
      cat("An error has occurred:\n")
      cat(e$message, "\n")
      # You can log the error or perform other actions
    }, warning = function(w) {
      # Code to run when there's a warning
      cat("Warning:\n")
      cat(w$message, "\n")
    }, finally = {
      validationRes<-list()
      for (i in 1:nrow(ligand_df)){
        row <- ligand_df[i,]
        upstreamCytokine<-row$ligand %>% as.character()
        receiver <- row$target_cell
        expressed_genes_receiver <- get_expressed_genes(receiver, ucint, pct = 0.10)
        background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
        ligands_of_interest <-c(upstreamCytokine,upstreamCytokine)
        gene_predictions_top20_list <- seq(5) %>%
          lapply(assess_rf_class_probabilities,
                 folds = 5, geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
                 ligands_oi = ligands_of_interest, ligand_target_matrix = ligand_target_matrix
          )

        # get top predicted genes
        top_predicted_genes <- seq(length(gene_predictions_top20_list)) %>%
          lapply(get_top_predicted_genes, gene_predictions_top20_list) %>%
          reduce(full_join, by = c("gene", "true_target"))
        top_predicted_cytokine_genes <- top_predicted_genes %>%
          filter(true_target) %>%
          dplyr::filter(gene %in% geneset_oi) %>% drop_na()
        top_predicted_cytokine_genes$upstrm<-upstreamCytokine
        top_predicted_cytokine_genes$target_cell <- receiver
        top_predicted_cytokine_genes<-top_predicted_cytokine_genes%>% dplyr::select('upstrm','gene', 'target_cell')
        print(top_predicted_cytokine_genes)
        validationRes[[paste0(upstreamCytokine,':',row$target_cell)]]<-as_tibble(top_predicted_cytokine_genes)
      }
      saveRDS(validationRes, paste0(output_file))
      cat("Finished executing the script.\n")
    }
  )
}

# UC
rf_validation_CKL(rds_file = 'uc_inflamed_naive.rds',
                  tsv_file = 'uc_inflamed_scibd_p01_weights_naive_new_model_input.tsv',
                  output_file = 'uc_inflamed_naive_validation.rds')

rf_validation_CKL(rds_file = 'uc_inflamed_no_naive.rds',
                  tsv_file = 'uc_inflamed_scibd_p01_weights_no_naive_new_model_input.tsv',
                  output_file = 'uc_inflamed_no_naive_validation.rds')

rf_validation_CKL(rds_file = 'uc_noninflamed_naive.rds',
                  tsv_file = 'uc_noninflamed_scibd_p01_weights_naive_new_model_input.tsv',
                  output_file = 'uc_noninflamed_naive_validation.rds')

rf_validation_CKL(rds_file = 'uc_noninflamed_no_naive.rds',
                  tsv_file = 'uc_noninflamed_scibd_p01_weights_no_naive_new_model_input.tsv',
                  output_file = 'uc_noninflamed_no_naive_validation.rds')


# Healthy
rf_validation_CKL(rds_file = 'scibd_healthy.rds',
                  tsv_file = 'healthy_scibd_p01_weights_new_model_input.tsv',
                  output_file = 'healthy_validation.rds')