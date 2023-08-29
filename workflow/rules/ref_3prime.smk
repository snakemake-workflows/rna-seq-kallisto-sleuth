if is_3prime_experiment:

    rule cds_polyA_T_removal:
        input:
            ref_fasta="resources/transcriptome.cdna.fasta",
        output:
            "resources/transcriptome.cdna.without_poly_a.fasta",
        log:
            "logs/kallisto_cds/cds_polyA_T_removal.log",
        conda:
            "../envs/biopython.yaml"
        script:
            "../scripts/remove_poly_tails.py"

    rule get_canonical_ids:
        output:
            "resources/canonical_ids.bed",
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
            fasta="resources/transcriptome.cdna.without_poly_a.fasta",
            canonical_ids="resources/canonical_ids.bed",
        output:
            "resources/transcriptome.cdna.without_poly_a.canonical.fasta",
        log:
            "logs/get_canonical_transcripts/get_canonical_transcripts.log",
        conda:
            "../envs/bedtools.yaml"
        shell:
            # TODO: make this use a BED file, possibly bedtools getfasta
            "bedtools getfasta "
            "  -fullHeader "
            "  -fi {input.fasta} "
            "  -bed {input.canonical_ids} "
            "  -fo {output} "
            "2> {log}"
