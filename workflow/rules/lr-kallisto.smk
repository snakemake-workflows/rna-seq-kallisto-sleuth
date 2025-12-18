# Example workflow
# https://github.com/pachterlab/LSRRSRLFKOTWMWMP_2024/blob/main/lr_kallisto_example.ipynb
# Paper
# https://pubmed.ncbi.nlm.nih.gov/39071335/
# Issue
# https://github.com/pachterlab/kallisto/issues/456

# --long options in the manual https://pachterlab.github.io/kallisto/manual -> -x bulk —long
# kallisto index will need run with -k 63.
# What would be the steps in order to generate the transcript-compatibility-counts-file required by quant-tcc -> The necessary files can all be generated from bustools count in bustools." So the commands would be kallisto bus, bustools sort, bustools correct (optional), bustools count, kallisto quant-tcc


#
# gffread -F -w GCA_000001405.15_GRCh38_no_alt_analysis_set_gencode_v45.fasta \
#    -g GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#    gencode.v45.annotation.gtf

# kallisto index -k 63 -t 10 -i gencode_v45 GCA_000001405.15_GRCh38_no_alt_analysis_set_gencode_v45.fasta

# Running lr-kallisto:

# kallisto bus -t 8 --long --threshold 0.8 -x bulk -i gencode_v45 \
#   -o kallisto_out fullLength.and.rescued.fastq

# bustools sort -t 8 kallisto_out/output.bus \
#  -o kallisto_out/sorted.bus

# bustools count kallisto_out/sorted.bus \
#  -t kallisto_out/transcripts.txt \
#  -e kallisto_out/matrix.ec \
#  -g kallisto_out/gencode_v45_tx2g.tsv \
#  -o kallisto_out/count --cm -m

# kallisto quant-tcc -t 8 \
# 	--long -p ONT -f kallisto_out/flens.txt \
# 	-i kallisto_index/gencode_v45 \
# 	-e kallisto_out/count.ec.txt \
# 	-o kallisto_out/quant-tcc \
# 	--matrix-to-files \
# 	kallisto_out/count.mtx


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
    log:
        "logs/kallisto_long_cdna/bus/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    threads: 5
    conda:
        "../envs/kallisto.yaml"
    shell:
        # Why -x? This is not single-cell?
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
    output:
        count_prefix="results/kallisto_long_cdna/{sample}-{unit}/count",
        matrix_ec="results/kallisto_long_cdna/{sample}-{unit}/matrix.ec",
        tx2g="results/kallisto_long_cdna/{sample}-{unit}/tx2g.tsv",
    log:
        "logs/kallisto_long_cdna/bustools_count/{sample}-{unit}.log",
    threads: 5
    conda:
        "../envs/bustools.yaml"
    shell:
        """
        bustools count {input.sorted_bus_file} \
        -t {input.transcripts} \
        -e {output.matrix_ec} \
        -g {output.tx2g} \
        -o {output.count_prefix} --cm -m 2> {log}
        """


# Was ist ${output}/count.mtx
rule kallisto_long_quant_tcc:
    input:
        # I have no equivalence classes file so I hope they can be taken from the index?
        # count_ec="results/kallisto_long_cdna/{sample}-{unit}/count.ec.txt",
        # -e {input.count_ec} \
        index="results/kallisto_long_cdna/transcripts.cdna.long.idx",
        # I have no fragment-file, so I do not perform normalization? Or should I use the fragement-length from the config?
        # flens="results/kallisto_long_cdna/{sample}-{unit}/flens.txt",
        # -f {input.flens}
    output:
        quant_folder=directory("results/kallisto_long_cdna/{sample}-{unit}/quant-tcc"),
    log:
        "logs/kallisto_long_cdna/quant_tcc/{sample}-{unit}.log",
    threads: 5
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto quant-tcc -t {threads} --long --platform ONT  \
        -i {input.index} \
        -o $(dirname {output.quant_folder}) \
        --matrix-to-files 2> {log}
        """
