rule get_transcriptome:
    output:
        "resources/transcriptome.{type}.fasta",
    log:
        "logs/get-transcriptome/{type}.log",
    params:
        species=config["resources"]["ref"]["species"],
        datatype="{type}",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
    wildcard_constraints:
        type="cdna|cds|ncrna",
    cache: True
    wrapper:
        "v1.7.1/bio/reference/ensembl-sequence"

if config["experiment"]["is-3-prime-rna-seq"]:

    rule cds_polyA_T_removal:
        input:
            ref_fasta="resources/transcriptome.cdna.fasta"
        output:
            "resources/transcriptome_clean.cdna.fasta",
        log:
            "results/logs/kallisto_cds/cds_polyA_T_removal.log",
        conda:
            "../envs/r-fasta.yaml"
        script:
            "../scripts/remove_poly_tails.py"

    rule get_3prime_seqs:
        input:
            read_length="results/stats/max-read-length.json",
            ref_fasta="resources/transcriptome_clean.cdna.fasta",
        output:
            "resources/transcriptome.3prime.fasta",
        params:
            release=config["resources"]["ref"]["release"],
        conda:
            "../envs/r-fasta.yaml"
        script:
            "../scripts/get_3prime-seqs.py"
    
    rule get_canonical_ids:
        output:
            "resources/canonical_ids.csv",
        log:
           "logs/filter_canonical/get_canonical_ids.log",
        params:
            release=config["resources"]["ref"]["release"],
        conda:
            "../envs/get_canonical_ids.yaml"
        script:
            "../scripts/get_canonical_ids.R"
    
    rule get_canonical_transcripts:
        input:
            fasta="resources/transcriptome.3prime.fasta",
            canonical_ids="resources/canonical_ids.csv",
        output:
            "resources/transcriptome_canonical.3prime.fasta",
        conda:
            "../envs/get_canonical_ids.yaml"
        shell:
            """bioawk -cfastx 'BEGIN{{while((getline k <"{input.canonical_ids}")>0)i[k]=1}}{{if(i[$name])print ">"$name"\\n"$seq}}' {input.fasta} > {output}"""

rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["resources"]["ref"]["species"],
        release=config["resources"]["ref"]["release"],
        build=config["resources"]["ref"]["build"],
        fmt="gtf",
    log:
        "logs/get-annotation.log",
    cache: True
    wrapper:
        "0.80.1/bio/reference/ensembl-annotation"


rule get_transcript_info:
    output:
        "resources/transcript-info.rds",
    params:
        species=get_bioc_species_name(),
        version=config["resources"]["ref"]["release"],
    log:
        "logs/get_transcript_info.log",
    conda:
        "../envs/biomart.yaml"
    cache: True
    script:
        "../scripts/get-transcript-info.R"


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
        "logs/convert-pfam.log",
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
    log:
        "logs/calculate-cpat-hexamers.log",
    conda:
        "../envs/cpat.yaml"
    cache: True
    shell:
        "make_hexamer_tab.py --cod={input.cds} --noncod={input.ncrna} > {output} 2> {log}"


rule calculate_cpat_logit_model:
    input:
        hexamers="resources/cpat.hexamers.tsv",
        cds="resources/transcriptome.cds.fasta",
        ncrna="resources/transcriptome.ncrna.fasta",
    output:
        "resources/cpat.logit.RData",
    params:
        prefix=lambda _, output: output[0][:-12],
    log:
        "logs/calculate-cpat-logit-model.log",
    conda:
        "../envs/cpat.yaml"
    cache: True
    shell:
        "make_logitModel.py --hex={input.hexamers} --cgene={input.cds} "
        "--ngene={input.ncrna} -o {params.prefix} 2> {log}"
