Code to reproduce scRNA-seq analysis in Fanti et al., 2022

The analysis performed is structured as follows (arrows indicate dependencies)

![](docs/dependencies.jpg)

The matching R environment is stored in `renv.lock` and can be restored using the `renv` package.

Raw and processed data are [available on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193322) and will be made publicly available upon publication. To reproduce the analysis, processed data should be placed in `data/day_<1/3>`.

Please make sure that other folders used throughout the analysis exist before running the scripts.
