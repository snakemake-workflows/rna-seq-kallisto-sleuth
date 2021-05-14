rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["resources"]["ref"]["species"],
        datatype="cdna",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"
