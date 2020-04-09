ensembl = config["resources"]["ref"]["ensembl"]

rule get_transcriptome:
    output:
        "results/refs/transcriptome.{type}.fasta"
    params:
        species=ensembl["species"],
        datatype="{type}",
        build=ensembl["build"],
        release=ensembl["release"]
    log:
        "logs/get-genome/{type}.log"
    wildcard_constraints:
        type="cdna|cds|ncrna"
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "results/refs/genome.gtf"
    params:
        species=ensembl["species"],
        release=ensembl["release"],
        build=ensembl["build"],
        fmt="gtf",
        flavor="chr_patch_hapl_scaff"
    wrapper:
        "0.50.0/bio/reference/ensembl-annotation"


rule get_pfam:
    output:
        r"results/refs/pfam/Pfam-A.{ext,(hmm|hmm\.dat)}"
    params:
        release=config["resources"]["ref"]["pfam"]
    log:
        "logs/get_pfam.{ext}.log"
    shell:
        "(curl -L ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{params.release}/Pfam-A.{wildcards.ext}.gz | gzip -d > {output}) 2> {log}"


rule convert_pfam:
    input:
        "results/refs/pfam/Pfam-A.hmm"
    output:
        multiext("results/refs/pfam/Pfam-A.hmm", ".h3m", ".h3i", ".h3f", ".h3p")
    log:
        "logs/convert_pfam.log"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmpress {input} > {log} 2>&1"


rule calculate_cpat_hexamers:
    input:
        cds="results/refs/transcriptome.cds.fasta",
        ncrna="results/refs/transcriptome.ncrna.fasta"
    output:
        "results/refs/cpat.hexamers.tsv"
    conda:
        "../envs/cpat.yaml"
    shell:
        "make_hexamer_tab.py --cod={input.cds} --noncod={input.ncrna} > {output}"


rule calculate_cpat_logit_model:
    input:
        hexamers="results/refs/cpat.hexamers.tsv",
        cds="results/refs/transcriptome.cds.fasta",
        ncrna="results/refs/transcriptome.ncrna.fasta"
    output:
        "results/refs/cpat.logit.RData"
    params:
        prefix=lambda _, output: output[0][:-12]
    conda:
        "../envs/cpat.yaml"
    shell:
        "make_logitModel.py --hex={input.hexamers} --cgene={input.cds} "
        "--ngene={input.ncrna} -o {params.prefix}"

