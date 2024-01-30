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


rule get_main_transcripts_fasta:
    input:
        fasta="resources/transcriptome.cdna.without_poly_a.fasta",
        main_transcripts_per_gene="resources/transcripts_annotation.main_transcript_strand_length.tsv",
    output:
        "resources/transcriptome.cdna.without_poly_a.main_transcript.fasta",
    log:
        "logs/get_canonical_transcripts/get_canonical_transcripts.log",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/get_main_transcripts_per_gene_fasta.py"
