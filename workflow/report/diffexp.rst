**Differentially expressed transcripts** using the model ``{{ snakemake.params.model["full"] }}``, computed with sleuth.
The columns ``b_*`` and ``b_*_se`` display the effect size :math:`\beta` and the corresponding standard error for every covariate in the model. 
They are analog to log2 fold changes. 
Note that non-binary covariates are binarized by sleuth, leading to multiple columns ``b_*1``, ``b_*2``, etc. 
In addition, a signed version of the pi-value score (as proposed by `Xiao et al. 2014 <https://dx.doi.org/10.1093/bioinformatics/btr671>`_) is shown. 
The sign reflects the sign of the effect (i.e. positive for upregulation, negative for downregulation).
