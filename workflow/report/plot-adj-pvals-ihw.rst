For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model }}``, raw p-values are plotted against adjusted p-values from IHW calculation for groups based on the covariate mean_obs.
Ideally, the curves show a similar trend, but with different degrees of adjustment of the p-value in the different groups.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#raw-versus-adjusted-p-values>`_ of the IHW package.
