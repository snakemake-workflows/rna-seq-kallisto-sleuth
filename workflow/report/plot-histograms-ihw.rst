For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model }}`` are shown a histogram of p-values overall and histograms for each group in IHW-grouping based on the covariate mean_obs.
Ideally, the histograms should show a uniform distribution at large p-values and any peak of p-values should be towards zero.
Only such a distribution ensures that IHW will control the false discovery rate.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#stratified-p-value-histograms>`_ of IHW package.
