if is_3prime_experiment:

    rule cds_polyA_T_removal:
        input:
            ref_fasta="resources/transcriptome.cdna.fasta",
        output:
            "resources/transcriptome_clean.cdna.fasta",
        log:
            "results/logs/kallisto_cds/cds_polyA_T_removal.log",
        conda:
            "../envs/biopython.yaml"
        script:
            "../scripts/remove_poly_tails.py"

    rule remove_strand_info_from_transcript_header:
        input:
            ref_fasta="resources/transcriptome_clean.cdna.fasta",
        output:
            "resources/transcriptome.3prime.fasta",
        params:
            release=config["resources"]["ref"]["release"],
        log:
            "logs/remove_strand_info_from_transcript_header/remove_strand_info_from_transcript_header.log",
        conda:
            "../envs/biopython.yaml"
        script:
            "../scripts/remove_strand_info.py"

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
            "resources/transcriptome_clean.3prime.fasta",
        log:
            "logs/get_canonical_transcripts/get_canonical_transcripts.log",
        conda:
            "../envs/get_canonical_ids.yaml"
        shell:
            """bioawk -cfastx \
            'BEGIN{{while((getline k <"{input.canonical_ids}")>0)i[k]=1}} \
            {{if(i[$name])print ">"$name"\\n"$seq}}' \
            {input.fasta} > {output}"""
