Diagnostic histograms for the mean of observations as choosen covariate in independent hypothesis weighting (IHW) calculation.
For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model["full"] }}`` there are displayed a histogram of p-values overall and histograms for each group in IHW-grouping.
Ideally the histograms should show a uniform distribution at large p-values. In this case IHW will control the false discovery rate.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#stratified-p-value-histograms>`_ of IHW package.
