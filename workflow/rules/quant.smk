rule kallisto_index:
    input:
        config["resources"]["ref"]["transcriptome"],
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
    shell:
        "kallisto quant -i {input.idx} -o {output} "
        "{params.extra} {input.fq} 2> {log}"
