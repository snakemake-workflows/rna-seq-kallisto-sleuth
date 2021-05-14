rule get_transcriptome:
    output:
        "resources/transcriptome.fasta",
    log:
        "logs/get-transcriptome.log",
    params:
        species=config["resources"]["ref"]["species"],
        datatype="cdna",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"
