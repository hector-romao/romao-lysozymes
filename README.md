[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code and data to reproduce the PGLS analyses in: "Duplication, Loss, and Feeding Ecology of Lysozyme Gene Copy Number in Insects"

---

## Repository structure
```
.
├──Supplementary_material #Supplementary figures in high quality
├──PGLS_analysis
├── data/                        # Input data files
├──── c_type_input_data.csv; i_type_input_data.csv: Input data with the regimes for each species used to the rescale the tree branches
├──── pgls_c_type_input_data_adults.csv; pgls_c_type_input_data_larave.csv; pgls_i_type_input_data_adults.csv; pgls_i_type_input_data_larave.csv: Input data with the gene count and diet for each species analysed
├──── ultrametric.nwk: Ultramerized tree
├── pgls_analysis.R              # Phylogenetic generalized least squares (PGLS) analysis
├── rescale_tree_C-type.R        # Phylogenetic tree rescaling – C-type
├── rescale_tree_I-type.R        # Phylogenetic tree rescaling – I-type
└── README.md
```

## Requirements

- **R** (≥ 4.0.0) — [Download](https://cran.r-project.org/)
- R packages: OUwie, phytools, ape and nlme

> Adjust the package list to match what your scripts actually use.

## How to reproduce the analyses

Run the scripts in the following order:

### 1. Tree rescaling
```r
Rscript rescale_tree_C-type.R
Rscript rescale_tree_I-type.R
```

These scripts rescale the phylogenetic tree for C-type and I-type analyses,
and save the rescaled trees to the working directory.

### 2. PGLS analysis
```r
Rscript pgls_analysis.R
```

Runs the PGLS regressions and outputs summary tables and figures.

## Data

The `data/` folder contains all input files necessary to run the scripts.
No additional downloads are required.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
