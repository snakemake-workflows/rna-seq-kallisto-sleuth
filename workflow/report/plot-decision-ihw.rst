Diagnostic plot for decision boundary of choosen covariate (mean of observations) after independent hypothesis weighting (IHW) calculation.
For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model["full"] }}`` there are displayed ranks of counts against unweighted p-values.
The coloured lines show the fold determined by IHW calculation and represent the decision boundary.
Ideally decisions are less tolerant at low covariate counts and require a higher p-value.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#decision-boundary>`_ of the IHW package.
