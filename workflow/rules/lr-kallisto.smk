# Example workflow from the paper (for single cell)
# https://github.com/pachterlab/LSRRSRLFKOTWMWMP_2024/blob/main/lr_kallisto_example.ipynb
# Paper
# https://pubmed.ncbi.nlm.nih.gov/39071335/
# Issue (for bulk)
# https://github.com/pachterlab/kallisto/issues/456


rule kallisto_long_index:
    input:
        fasta=(
            "resources/transcriptome.cdna.without_poly_a.fasta"
            if is_3prime_experiment
            else "resources/transcriptome.cdna.fasta"
        ),
    output:
        index="results/kallisto_long_cdna/transcripts.cdna.long.idx",
    log:
        "logs/kallisto_long_cdna/index.cdna.log",
    threads: 10
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto index -k 63 -t {threads} -i {output.index} {input.fasta} > {log} 2>&1
        """


rule kallisto_long_bus:
    input:
        fastq=kallisto_quant_input,
        index="results/kallisto_long_cdna/transcripts.cdna.long.idx",
    output:
        bus_file="results/kallisto_long_cdna/{sample}-{unit}/output.bus",
        transcript="results/kallisto_long_cdna/{sample}-{unit}/transcripts.txt",
        flens="results/kallisto_long_cdna/{sample}-{unit}/flens.txt",
        matrix="results/kallisto_long_cdna/{sample}-{unit}/matrix.ec",
    log:
        "logs/kallisto_long_cdna/bus/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto bus --threshold 0.8 --long -x bulk -t {threads} -i {input.index} \
        -o $(dirname {output.bus_file}) {input.fastq} 2> {log}
        """


rule bustools_sort:
    input:
        bus_file="results/kallisto_long_cdna/{sample}-{unit}/output.bus",
    output:
        sorted_bus_file="results/kallisto_long_cdna/{sample}-{unit}/sorted.bus",
    log:
        "logs/kallisto_long_cdna/bustools_sort/{sample}-{unit}.log",
    threads: 5
    conda:
        "../envs/bustools.yaml"
    shell:
        """
        bustools sort -t {threads} {input.bus_file} \
        -o {output.sorted_bus_file} 2> {log}
        """


rule bustools_count:
    input:
        sorted_bus_file="results/kallisto_long_cdna/{sample}-{unit}/sorted.bus",
        transcripts="results/kallisto_long_cdna/{sample}-{unit}/transcripts.txt",
        matrix_ec="results/kallisto_long_cdna/{sample}-{unit}/matrix.ec",
        transcript_info="resources/transcripts_annotation.results.tsv",
    output:
        # We only need this file for computation with default-storage-provider
        tmp_file=temp("results/kallisto_long_cdna/{sample}-{unit}/count"),
        count_prefix="results/kallisto_long_cdna/{sample}-{unit}/count.ec.txt",
        count_mtx="results/kallisto_long_cdna/{sample}-{unit}/count.mtx",
    log:
        "logs/kallisto_long_cdna/bustools_count/{sample}-{unit}.log",
    threads: 5
    conda:
        "../envs/bustools.yaml"
    params:
        prefix="results/kallisto_long_cdna/{sample}-{unit}/count",
    shell:
        """
        bustools count {input.sorted_bus_file} \
        -t {input.transcripts} \
        -e {input.matrix_ec} \
        -g {input.transcript_info} \
        -o {output.tmp_file} --cm -m 2> {log}
        touch {output.tmp_file}
        """


rule kallisto_long_quant_tcc:
    input:
        count_ec="results/kallisto_long_cdna/{sample}-{unit}/count.ec.txt",
        count_mtx="results/kallisto_long_cdna/{sample}-{unit}/count.mtx",
        index="results/kallisto_long_cdna/transcripts.cdna.long.idx",
        flens="results/kallisto_long_cdna/{sample}-{unit}/flens.txt",
    output:
        quant_folder=directory("results/kallisto_long_cdna/{sample}-{unit}/quant-tcc"),
    log:
        "logs/kallisto_long_cdna/quant_tcc/{sample}-{unit}.log",
    threads: 5
    conda:
        "../envs/kallisto.yaml"
    params:
        platform=config["sequencing_platform"],
    shell:
        """
        kallisto quant-tcc -t {threads} -b 3 --long --platform {params.platform}  \
        -e {input.count_ec} \
        -f {input.flens} \
        -i {input.index} \
        -o {output.quant_folder} \
        --matrix-to-files {input.count_mtx} 2> {log} 
        cp {output.quant_folder}/abundance_1.h5 {output.quant_folder}/abundance.h5
        """