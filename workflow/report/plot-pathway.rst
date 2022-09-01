**Pathway enrichment** performed with SPIA, using the model ``{{ snakemake.config["diffexp"]["models"][snakemake.wildcards.model]["full"] }}``.

The plot is filtered based on ``pGFdr`` <0.05   On the x-axis contains the ``tA`` is the observed total perturbation accumulation in the pathway and y-axis contains the pathways. Bars are coloured based on the ``Status`` gives the direction in which the pathway is perturbed (``activated`` or ``inhibited``) (also see the `SPIA docs <https://rdrr.io/bioc/SPIA/man/spia.html>`_).

