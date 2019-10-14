After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_, transcripts were quantified with `Kallisto <https://pachterlab.github.io/kallisto/>`_.
Integrated normalization and differential expression analysis was conducted with `Sleuth <https://pachterlab.github.io/sleuth>`_ following standard procedure as outlined in the manual.
For sample metadata, see {{ snakemake.config["samples"] }}_.
