Histogram of differential expression p-values computed with sleuth using the model ``{{ snakemake.wildcards.model }}``, aggregated to the level of {{ snakemake.wildcards.level }}.

For the false discovery rate (FDR) control to be valid, the p-value distribution must be `super-uniform <https://stats.stackexchange.com/q/419005/254369>`_.
Roughly speaking, this means that the distribution should fall monotonically from zero towards one and a peak at exactly one is the only tolerable increase.
For a detailed definition have a look at `this cross-validated discussion of the definition <https://stats.stackexchange.com/q/419005/254369>`_ and the following article as an example where this constraint is formulated as ``stochastically larger than uniform``:
Benjamini, Y., Krieger, A.M., and Yekutieli, D. (2006). Adaptive linear step-up procedures that control the false discovery rate. Biometrika 93, 491–507 . doi:10.1093/biomet/93.3.491.
