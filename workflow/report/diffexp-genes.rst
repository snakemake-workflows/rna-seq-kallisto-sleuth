**Differentially expressed genes** using the model ``{{ snakemake.params.model["full"] }}``, computed with sleuth by aggregating transcript p-values.
In addition, a signed version of the pi-value score (as proposed by `Xiao et al. 2014 <https://dx.doi.org/10.1093/bioinformatics/btr671>`_) is shown.
The sign reflects the sign of the effect (i.e. positive for upregulation, negative for downregulation).
