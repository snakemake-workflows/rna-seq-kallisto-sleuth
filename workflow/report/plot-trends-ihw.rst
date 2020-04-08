Diagnostic plot of general trend of estimated weights of mean of observations (choosen covariate) after independent hypothesis weighting (IHW) calculation.
For each ``{{ snakemake.wildcards.level }}`` of ``{{ snakemake.wildcards.model["full"] }}`` there are displayed trends coloured by different folds after IHW calculation.
Ideally the curves should show a similar trend with a priorisation in weighting at high means of counts.
For more information see `documentation <https://www.bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#estimated-weights>`_ of the IHW package.
