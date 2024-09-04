# Function: iterateCells
# -----------------------
# A function that takes sender cell types, a receiver cell type, and a dataset as inputs. 
# The function filters for cytokine ligands based on test results and predicts their activities using 
# NicheNet's ligand-target matrix. The function returns a list of sender cell types, 
# the receiver cell type, and the cytokine ligands that are significant based on either 
# AUROC or Pearson correlation coefficient.
#
# Parameters:
# -----------
# sndr: A character vector of sender cell types
# rcvr: A character vector of a single receiver cell type
# dataset: A data frame that contains gene expression data of the sender and receiver cells
#
# Returns:
# --------
# A list that contains the following components:
# 1. A character vector of sender cell types
# 2. A character vector of a single receiver cell type
# 3. A character vector of significant cytokine ligands based on either AUROC or Pearson correlation coefficient.
#
# Note: If no significant cytokine ligands are found, the function prints a message to the console and returns NULL.

iterateCells <- function(sndr,rcvr,dataset) {
  
  tryCatch(
    expr = {
      ## sender
      # Get unique sender cell types
      sender_celltypes <- sndr
      
      # Get the expressed genes of every sender cell type separately
      list_expressed_genes_sender <- sender_celltypes %>%
        unique() %>%
        lapply(get_expressed_genes, dataset, 0.10)
      
      # Combine the expressed genes of all the sender cell types
      expressed_genes_sender <- list_expressed_genes_sender %>%
        unlist() %>%
        unique()
      
      ## receiver
      # Get expressed genes of the receiver cell type
      # using the 10% cutoff for 10X data as recommended by the NicheNet authors 
      # https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md
      receiver <- rcvr
      expressed_genes_receiver <- get_expressed_genes(receiver, dataset, pct = 0.10)
      
      # Get background expressed genes that are included in NicheNet's ligand-target matrix
      background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
      
      ligands = lr_network %>% pull(ligand) %>% unique()%>% .[. %in% colnames(ligand_target_matrix)] %>% .[. %in% geneset_oi]
      expressed_ligands = intersect(ligands,expressed_genes_sender)
      receptors = lr_network %>% pull(receptor) %>% unique()
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      lr_network_expressed = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) 
      potential_ligands = lr_network_expressed %>% pull(ligand) %>% unique()
      # Predict the activities of the cytokine ligands using NicheNet's ligand-target matrix
      nichenet_activities <- predict_ligand_activities(
        geneset = geneset_oi,
        background_expressed_genes = background_expressed_genes,
        ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands
      ) 
      
      # Filter for significant cytokine ligands based on either AUROC or Pearson correlation coefficient
      best_upstream_ligands<- nichenet_activities%>% dplyr::filter(pearson >= 0.1) %>% pull(test_ligand)
      active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
      active_ligand_target_links_df$source_cell <- paste0(i)
      active_ligand_target_links_df$target_cell <- paste0(j)
      # If significant cytokine ligands are found, return 
      # the list of sender cell types, 
      # receiver cell type, 
      # and significant cytokine ligands
      if (length(best_upstream_ligands)> 0) {
        #print(paste0(sndr," - ",rcvr))
        #return(nichenet_activities)
        return(active_ligand_target_links_df)
      } else if (length(ligands_of_interest) == 0) {
        #print(paste0("No significant cytokine ligands found for ",sndr," - ",rcvr))
      }
    },
    error = function(e){ 
      print(e)
    },
    warning = function(w){
      print(w)
    },
    finally = {
      # (Optional)
      # Do this at the end before quitting the tryCatch structure...
    }
  )
}