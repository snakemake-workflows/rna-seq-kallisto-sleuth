rule cutadapt_pe:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        "{}".format(config["params"]["cutadapt-pe"]),
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "0.31.1/bio/cutadapt/pe"


if is_3prime_experiment:

    # Here rule cutadapt1 checks and remove poly-A tails and sequence qualilty score <20.
    # reasoning behind parameters:
    #   * `-m 20`:
    #   * Discards any read under 20bp by -m 20
    #   * `-O 20 -a ""polyA=A{20}"" -a ""QUALITY=G{20}"" -n 2`:
    #   * Removes reads having polyA (18bp) or polyG (18bp) and repeat the process twice by -n 2
    # In short, it removes reads that are predominantly polyA or polyG.
    rule cutadapt1:
        input:
            get_fastqs,
        output:
            fastq="results/trimmed/{sample}-{unit}.1.1.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.1.1.qc.txt",
        params:
            extra=lambda w: str(
                "-m 20 -O 20 -a " "polyA=A{20}" " -a " "QUALITY=G{20}" " -n 2"
            ),
        log:
            "results/logs/cutadapt/{sample}-{unit}.1.1.log",
        wrapper:
            "v1.14.1/bio/cutadapt/se"

    # Here rule cutadapt2 checks if reads are obtained from nextseq and removes adapters specific to nextseq
    # reasoning behind parameters:
    #   * `-m 20`:
    #   * Discards any read under 20bp by -m 20
    #   * `--nextseq-trim=10`
    #   * Quality score<20 is trimmed, based on NextSeq/NovaSeq 2-color SBS system
    #   * `-a A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000`
    #   * Trims the polyA-stretch + adapter at the 3' end with minimum overlap 3bp, maximum error of 1bp every 10bp
        rule cutadapt2:
            input:
                fastq="results/trimmed/{sample}-{unit}.1.1.fastq.gz",
            output:
                fastq="results/trimmed/{sample}-{unit}.2.1.fastq.gz",
                qc="results/trimmed/{sample}-{unit}.2.1.qc.txt",
            params:
                adapters=lambda w: str(
                    '-m 20 -O 3 --nextseq-trim=10 -a "'
                    'A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000"'
                    ""
                ),
            log:
                "results/logs/cutadapt/{sample}-{unit}.2.1.log",
            wrapper:
                "v1.14.1/bio/cutadapt/se"

    # Here rule cutadapt3 removes final set of adapters from 5'end and discards reads based on trimmed read length
    # reasoning behind parameters:
    #   * `-m 20`:
    #   * Discards any read under 20bp by -m 20
    #   * `-O 20 -g r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20 --discard-trimmed"`
    #   * Trim the 5' end of adapter with miminum overlap of 20bp and discard its read by --discard-trimmed
        rule cutadapt3:
            input:
                fastq="results/trimmed/{sample}-{unit}.2.1.fastq.gz",
            output:
                fastq="results/trimmed/{sample}-{unit}.fastq.gz",
                qc="results/trimmed/{sample}-{unit}.qc.txt",
            params:
                adapters=lambda w: str(
                    '-m 20 -O 20 -g "'
                    'r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20"'
                    " --discard-trimmed"
                ),
            log:
                "results/logs/cutadapt/{sample}-{unit}.log",
            wrapper:
                "v1.14.1/bio/cutadapt/se"
else:

    rule cutadapt:
        input:
            get_fastqs,
        output:
            fastq="results/trimmed/{sample}-{unit}.fastq.gz",
            qc="results/trimmed/{sample}-{unit}.qc.txt",
        params:
            "{}".format(config["params"]["cutadapt-se"]),
        log:
            "results/logs/cutadapt/{sample}-{unit}.log",
        wrapper:
            "0.31.1/bio/cutadapt/se"


rule max_read_length:
    input:
        get_all_fastqs,
    output:
        "results/stats/max-read-length.json",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get-max-read-length.py"
