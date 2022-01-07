Code to reproduce scRNA-seq analysis in Fanti et al., 2022

The analysis performed is structured as follows (arrows indicate dependencies)

![](docs/dependencies.jpg)

The matching R environment is stored in `renv.lock` and can be restored using the `renv` package.

A GEO accession number will be added as soon as available. Processed data downloaded from GEO should be placed in `data/day_<1/3>`. Make sure that other folders used throughout the analysis exist before running the scripts.
