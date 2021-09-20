from pathlib import Path

ensembl = config["resources"]["ref"]["ensembl"]


rule get_transcriptome:
    output:
        "resources/transcriptome.{type}.fasta",
    params:
        species=ensembl["species"],
        datatype="{type}",
        build=ensembl["build"],
        release=ensembl["release"],
    log:
        "logs/get-genome/{type}.log",
    wildcard_constraints:
        type="cdna|cds|ncrna",
    cache: True
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=ensembl["species"],
        release=ensembl["release"],
        build=ensembl["build"],
        fmt="gtf",
        flavor="chr_patch_hapl_scaff",
    cache: True
    wrapper:
        "0.50.0/bio/reference/ensembl-annotation"


rule get_pfam:
    output:
        r"resources/pfam/Pfam-A.{ext,(hmm|hmm\.dat)}",
    params:
        release=config["resources"]["ref"]["pfam"],
    log:
        "logs/get_pfam.{ext}.log",
    shell:
        "(curl -L ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/"
        "Pfam{params.release}/Pfam-A.{wildcards.ext}.gz | "
        "gzip -d > {output}) 2> {log}"


rule convert_pfam:
    input:
        "resources/pfam/Pfam-A.hmm",
    output:
        multiext("resources/pfam/Pfam-A.hmm", ".h3m", ".h3i", ".h3f", ".h3p"),
    log:
        "logs/convert_pfam.log",
    conda:
        "../envs/hmmer.yaml"
    cache: True
    shell:
        "hmmpress {input} > {log} 2>&1"


rule calculate_cpat_hexamers:
    input:
        cds="resources/transcriptome.cds.fasta",
        ncrna="resources/transcriptome.ncrna.fasta",
    output:
        "resources/cpat.hexamers.tsv",
    conda:
        "../envs/cpat.yaml"
    cache: True
    shell:
        "make_hexamer_tab.py --cod={input.cds} --noncod={input.ncrna} > {output}"


rule calculate_cpat_logit_model:
    input:
        hexamers="resources/cpat.hexamers.tsv",
        cds="resources/transcriptome.cds.fasta",
        ncrna="resources/transcriptome.ncrna.fasta",
    output:
        "resources/cpat.logit.RData",
    params:
        prefix=lambda _, output: output[0][:-12],
    conda:
        "../envs/cpat.yaml"
    cache: True
    shell:
        "make_logitModel.py --hex={input.hexamers} --cgene={input.cds} "
        "--ngene={input.ncrna} -o {params.prefix}"


rule download_bioconductor_species_database:
    output:
        directory("resources/bioconductor/lib/R/library/{package}"),
    params:
        path=lambda wc, output: Path(output[0]).parents[3],
        version=config["resources"]["ref"]["species_db_version"],
    log:
        "logs/resources/bioconductor/{package}.log",
    shell:
        "conda create --yes --quiet -p {params.path} "
        "--channel bioconda bioconductor-{wildcards.package}={params.version}"
