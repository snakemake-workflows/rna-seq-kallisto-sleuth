For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model }}`` p-values are plotted against the ranks of counts of the covariate mean_obs.
Ideally the dispersion shows a trend, that describes the correlation of the covariate under the alternative hypothesis, for example low p-values are enriched at high covariate values.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#scatter-plots>`_ of the IHW package.
