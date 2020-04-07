Diagnostic plot about adjustment of p-values after independent hypothesis weighting (IHW) calculation for mean of observations as choosen covariate.
For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.params.model["full"] }}`` there are displayed raw p-values against adjusted p-values from IHW calculation.
Ideally, the curves show a similar trend, but with different degrees of adjustment of the p-value in the different groups.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#raw-versus-adjusted-p-values>`_ of the IHW package.
