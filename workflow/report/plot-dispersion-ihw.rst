Diagnostic scatter plot for the mean of observations as choosen covariate in independent hypothesis weighting (IHW) calculation.
For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model["full"] }}`` there are displayed the p-values vs. the ranks of counts.
Ideally the dispersion shows a trend, which describe the correlation of the covariate under the alternative hypothesis, for example low p-values are enriched at high covariate values.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#scatter-plots>`_ of the IHW package.
