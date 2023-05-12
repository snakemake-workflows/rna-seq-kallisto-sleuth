rule kallisto_index:
    input:
        fasta="resources/transcriptome_clean.cdna.fasta"
        if is_3prime_experiment
        else "resources/transcriptome.cdna.fasta",
    output:
        index="results/kallisto_cdna/transcripts.cdna.idx",
    log:
        "results/logs/kallisto_cdna/index.cdna.log",
    threads: 1
    wrapper:
        "v1.23.1/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=kallisto_quant_input,
        index="results/kallisto_cdna/transcripts.cdna.idx",
    output:
        kallisto_folder=directory("results/kallisto_cdna/{sample}-{unit}"),
    log:
        "results/logs/kallisto_cdna/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    wrapper:
        "v1.23.1/bio/kallisto/quant"
