rule fastp_se:
    input:
        sample=[ get_fastqs ],
    output:
        trimmed="results/trimmed/{sample}/{sample}-{unit}.fastq.gz",
        failed="results/trimmed/{sample}/{sample}-{unit}.failed.fastq.gz",
        html="results/trimmed/{sample}/{sample}-{unit}.html",
        json="results/trimmed/{sample}/{sample}-{unit}.json",
    log:
        "logs/trimmed/{sample}/{sample}-{unit}.log",
    params:
        adapters=lookup(within=units, query="sample == '{sample}' & unit =='{unit}'", cols="fastp_adapters", default=""),
        extra=lookup(within=units, query="sample == '{sample}' & unit =='{unit}'", cols="fastp_extra", default="--trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7 --length_required 33"),
    threads: 1
    wrapper:
        "v7.1.0/bio/fastp"


rule fastp_pe:
    input:
        sample=get_fastqs,
    output:
        trimmed=[
            "results/trimmed/{sample}/{sample}-{unit}.1.fastq.gz",
            "results/trimmed/{sample}/{sample}-{unit}.2.fastq.gz",
        ],
        # Unpaired reads separately
        unpaired1="results/trimmed/{sample}/{sample}-{unit}.unpaired.1.fastq.gz",
        unpaired2="results/trimmed/{sample}/{sample}-{unit}.unpaired.u2.fastq.gz",
        failed="results/trimmed/{sample}/{sample}-{unit}.failed.fastq.gz",
        html="results/trimmed/{sample}/{sample}-{unit}.html",
        json="results/trimmed/{sample}/{sample}-{unit}.json",
    log:
        "logs/trimmed/{sample}/{sample}-{unit}.log",
    params:
        adapters=lookup(within=units, query="sample == '{sample}' & unit =='{unit}'", cols="fastp_adapters", default="--detect_adapter_for_pe"),
        extra=lookup(within=units, query="sample == '{sample}' & unit =='{unit}'", cols="fastp_extra", default="--trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7 --length_required 33"),
    threads: 2
    wrapper:
        "v7.1.0/bio/fastp"


rule max_read_length:
    input:
        get_all_fastqs,
    output:
        "results/stats/max-read-length.json",
    log:
        "logs/max-read-length.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get-max-read-length.py"
