#Just replace the keys and values of datalist with your rds and tsv files. The output will be an edge list, stored in an RDS file, with source and target cytokines, collapsed mediating cell-cell interactions, n (number of cell-cell interactions) and the state(s if there are multiple).

# combine validated and original networks
library(tidyverse)


# list including file pairs to combine, match corresponding validated (rds) and original (tsv) networks
datalist<- list("validated_1.rds" = "1_cytokine_network.tsv",
                "validated_2.rds" = "2_cytokine_network.tsv"
                )
# list to save results in
results <- list()

# loop iterates over file pairs and saves the validated hits in a corresponding element of results
for (p in names(datalist)){
  statename <- paste0(datalist[[p]])  |> 
    str_replace(pattern = '_network.tsv', replacement = '') # you may need to change this, depending on the naming of your files
  orig <- read_tsv(paste0(datalist[[p]]))|> rename(source = ligand)
  valid <- readRDS(paste0(p)) |> bind_rows() |> rename(source = upstrm, target = gene)
  net <- left_join(valid, orig) |> drop_na()
  net_collapse<-net |> tidyr::unite(c('source_cell','target_cell'), col='CCC', sep=':') |> 
    group_by(source,target) |> 
    summarise(cells=paste(CCC, collapse=","), n= n())
  net_collapse$state <- paste0(statename)
  results[[paste0(statename)]] <- net_collapse
}

# join all files
outfile <- bind_rows(results)
# write them to file
saveRDS(outfile, '../combined_cytokine_network.rds')