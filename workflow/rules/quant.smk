rule kallisto_index:
    input:
        "resources/transcriptome.cdna.fasta",
    output:
        "results/kallisto/transcripts.idx",
    log:
        "results/logs/kallisto/index.log",
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"


rule kallisto_quant:
    input:
        fq=get_trimmed,
        idx="results/kallisto/transcripts.idx",
    output:
        directory("results/kallisto/{sample}-{unit}"),
    log:
        "results/logs/kallisto/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    # around 4 gb of memory usage with 1 thread. (hg38)
    # around 5.5gb peak memory usage with 8 threads. (hg38)
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        "kallisto quant -i {input.idx} -o {output} "
        "{params.extra} -t {threads} {input.fq} 2> {log}"
