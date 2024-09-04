# Cytokine Networks Reveal Pathogenic Mechanisms and Therapeutic Targets in Ulcerative Colitis

This repository contains the necessary scripts to replicate the analysis detailed in Olbei et al.:
- Olbei M., Hautefort I., Thomas P. J., Csabai L., Bohar B., Ibraheim H., Saifuddin A., Cozzetto D., Modos D., Powell N., Korcsmaros T.

The `codes` folder contains two subdirectories: `codes/ data_processing` refers to scripts handling the scRNA-Seq data and generation of the cytokine networks, while `codes/ figures` contains an Rmarkdown notebook detailing how to generate the figures shown in the manuscript.

## Data Processing
Brief description of main data processing scripts:
### `process_scibd.R`: 
- **Input**: Pre-processed single-cell RNA-Seq data downloaded from [scibd.cn](http://scibd.cn/), in `.rds` format. 
- **Task**: Separates the analysed sample groups (inflamed-treated, inflamed-naive, noninflamed-treated, healthy) into corresponding .rds files using metadata from the original publications included in scIBD. 
- **Output**: State-specific `.rds.` files

### `generate_models.R`: 
- **Input**: State-specific `.rds` files; nichenet prior models; geneset of interest (list of cytokines)
- **Task**: Generation of cytokine-cytokine interactions using nichenet for each state
- **Output**: Cytokine networks
- **Uses**: `iteratecells.R`: function containing the individual steps of the nichenet analysis, gets called by `generate_models.R`

### `nichenet_validation.R`: 
- **Input**: Cytokine networks
- **Task**: Nichenet's target prediction evaluation pipeline using a multi-ligand random forest model with cross validation. Target downstream cytokines from the gene set of interest are kept if they are well predicted in every cross-validation round
- **Output**: Validated cytokine networks
