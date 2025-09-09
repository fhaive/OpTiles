# OpTiles: An R Package for Adaptive Tiling and Methylation Variability Profiling

OpTiles, an R package that dynamically defines tiling windows based on the distribution of sequenced CpGs. OpTiles further quantifies intra-region methylation variability and integrates it with CpG density to generate a metric for region reliability and biological interpretability. The package also includes functions for annotating genomic regions and to quantify overlaps with regions of interest. While compatible with existing tools such as methylKit, OpTiles extends functionality by refining region definitions, supporting data-driven prioritization, and enhancing the interpretation of heterogeneous or complex datasets.

![Architecture Diagram](how_to_use/graphical_abstract_manual.png)

The package implements the methods and algorithms described in the following scientific publication:

> **G. Migliaccio** (2025). *OpTiles: An R Package for Adaptive Tiling and Methylation Variability Profiling*. bioRxiv BIORXIV/2025/674166, DOI: [].

If you use this package in your research, please cite the paper using the citation provided.

## Installation
R version >4.4.1
To install OpTiles, and download the packages

```
git clone https://github.com/fhaive/OpTiles.git
cd OpTiles

R
renv::restore()
devtools::install(".")

```

## Example script
The script `how_to_use/OpTiles_example_script.R` contains the steps to reproduce the analysis described in the manual.
All the input data used and the outputs produced are available on [Zenodo].
