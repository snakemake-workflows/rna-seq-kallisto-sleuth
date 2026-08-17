rule get_transcriptome:
    output:
        "resources/transcriptome.{type}.fasta",
    log:
        "logs/get-transcriptome/{type}.log",
    wildcard_constraints:
        type="cdna|cds|ncrna",
    cache: "omit-software"
    localrule: True
    params:
        species=config["resources"]["ref"]["species"],
        datatype="{type}",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
    wrapper:
        "v1.7.1/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    log:
        "logs/get-annotation.log",
    cache: "omit-software"
    localrule: True
    params:
        species=config["resources"]["ref"]["species"],
        release=config["resources"]["ref"]["release"],
        build=config["resources"]["ref"]["build"],
        flavor="chr_patch_hapl_scaff",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    wrapper:
        "v6.0.1/bio/reference/ensembl-annotation"


rule get_transcript_info:
    output:
        multiext(
            "resources/transcripts_annotation",
            ".results.rds",
            ".results.tsv",
            ".main_transcript_strand_length.tsv",
        ),
    log:
        "logs/get_transcript_info.log",
    cache: "omit-software"
    conda:
        "../envs/biomart.yaml"
    params:
        species=get_bioc_species_name(),
        version=config["resources"]["ref"]["release"],
        three_prime_activated=is_3prime_experiment,
    script:
        "../scripts/get-transcript-info.R"


rule get_pfam:
    output:
        r"resources/pfam/Pfam-A.{ext,(hmm|hmm\.dat)}",
    log:
        "logs/get_pfam.{ext}.log",
    cache: True
    localrule: True
    params:
        release=config["resources"]["ref"]["pfam"],
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
        "logs/convert-pfam.log",
    cache: True
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmpress {input} > {log} 2>&1"


rule calculate_cpat_hexamers:
    input:
        cds="resources/transcriptome.cds.fasta",
        ncrna="resources/transcriptome.ncrna.fasta",
    output:
        "resources/cpat.hexamers.tsv",
    log:
        "logs/calculate-cpat-hexamers.log",
    cache: True
    conda:
        "../envs/cpat.yaml"
    shell:
        "make_hexamer_tab.py --cod={input.cds} --noncod={input.ncrna} > {output} 2> {log}"


rule calculate_cpat_logit_model:
    input:
        hexamers="resources/cpat.hexamers.tsv",
        cds="resources/transcriptome.cds.fasta",
        ncrna="resources/transcriptome.ncrna.fasta",
    output:
        "resources/cpat.logit.RData",
    log:
        "logs/calculate-cpat-logit-model.log",
    cache: True
    conda:
        "../envs/cpat.yaml"
    params:
        prefix=lambda _, output: output[0][:-12],
    shell:
        "make_logitModel.py --hex={input.hexamers} --cgene={input.cds} "
        "--ngene={input.ncrna} -o {params.prefix} 2> {log}"


rule get_spia_db:
    input:
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        "resources/spia-db.{database}.rds",
    log:
        "logs/spia-db.{database}.log",
    cache: True
    retries: 3
    conda:
        enrichment_env
    params:
        bioc_species_pkg=bioc_species_pkg,
        species=get_bioc_species_name(),
    script:
        "../scripts/get-spia-db.R"


rule get_decon_references:
    output:
        rda="resources/celltype_references/{celltype_reference}.xCell2Ref.rda",
    log:
        "logs/xcell2/{celltype_reference}.get_decon_ref.log",
    localrule: True
    params:
        base_url="https://raw.githubusercontent.com/AlmogAngel/xCell2/master/data/",
    shell:
        """
        curl -fsSL {params.base_url}/{wildcards.celltype_reference}.xCell2Ref.rda -o {output.rda} >{log} 2>&1
        """
