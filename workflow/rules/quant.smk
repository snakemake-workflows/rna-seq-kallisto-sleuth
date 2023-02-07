rule kallisto_index:
    input:
        fasta="resources/transcriptome_clean.{type}.fasta"
        if is_3prime_experiment
        else "resources/transcriptome.cdna.fasta",
    output:
        index="results/kallisto_{type}/transcripts.{type}.idx",
    log:
        "results/logs/kallisto_{type}/index.{type}.log",
    wildcard_constraints:
        type="cdna|3prime",
    threads: 1
    wrapper:
        "v1.17.4/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=kallisto_quant_input,
        index="results/kallisto_{type}/transcripts.{type}.idx",
    output:
        kallisto_folder=directory("results/kallisto_{type}/{sample}-{unit}"),
    log:
        "results/logs/kallisto_{type}/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    wildcard_constraints:
        type="cdna|3prime",
    wrapper:
        "v1.17.4/bio/kallisto/quant"
