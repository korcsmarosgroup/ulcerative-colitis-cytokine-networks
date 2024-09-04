# processing scIBD data for the cytokine networks
library(Seurat)
library(tidyseurat)
library(tidyverse)
library(nichenetr)
library(OmnipathR)

# Loading data as downloaded from scIBD: http://scibd.cn/
adata <- readRDS('../data/scIBD.gex_matrix.rds')

# Selecting UC and healty data to help downstream processing
uc <- adata |> tidyseurat::filter(stage == 'adult') |> tidyseurat::filter(disease == 'UC_inflamed' | disease == 'UC_non_inflamed' | disease == 'Healthy')
saveRDS(uc,'../rds/uc_scibd.rds')
rm(adata)

# 7 inflamed naive patient vector, 3 from parikh, 2 from corridoni, 2 from kinchen
pat <- c('GSM3214201_A3', 'GSM3214204_B3', 'GSM3214207_C3', 'GSM4483695_S24','GSM4483696_S33','GSM3140595_UC1','GSM3140596_UC2')
# 3 noninflamed naive
ni_naive <- c("GSM3214202_A2", "GSM3214205_B2", "GSM3214208_C2")

ucint <- uc |> tidyseurat::filter(tissue == 'largeInt' | tissue == 'smallInt') 
rm(uc)
# =====
# UC inflamed
# UC inflamed naive
uci <- ucint |> tidyseurat::filter(disease == 'UC_inflamed')|> tidyseurat::filter(subject %in% pat)
saveRDS(uci, '../rds/uc_inflamed_naive.rds')
# UC inflamed (without naive patients)
uci_no_naive <- uci |> tidyseurat::filter(!subject %in% pat)
saveRDS(uci_no_naive, '../rds/uc_inflamed_no_naive.rds')

# =====
# UC non inflamed
ucn <- ucint |> tidyseurat::filter(disease == 'UC_non_inflamed')
rm(ucint)
# UC non-inflamed (without naive patients)
ucn_no_naive <- ucn |> tidyseurat::filter(!subject %in% ni_naive)
saveRDS(ucn_no_naive, '../rds/uc_noninflamed_no_naive.rds')

# UC non-inflamed naive (was not used in manuscript due to its small size)
#ucn_naive <- ucn |> tidyseurat::filter(subject %in% ni_naive)
#saveRDS(ucn_naive, '../rds/uc_noninflamed_naive.rds')

# =====
# healthy data
h <- ucint |> tidyseurat::filter(disease == 'Healthy')
saveRDS(h, '../rds/scibd_healthy.rds')