rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

