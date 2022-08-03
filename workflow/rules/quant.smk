rule kallisto_index:
    input:
        "resources/transcriptome.3prime.fasta"
        if config["experiment"]["is-3-prime-rna-seq"]
        else "resources/transcriptome.cdna.fasta",
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
        kallisto_folder=directory("results/kallisto/{sample}-{unit}"),
        bam_file="results/kallisto/{sample}-{unit}/pseudoalignments.bam",
    log:
        "results/logs/kallisto/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output.kallisto_folder} {params.extra} {input.fq} 2> {log}"


rule get_heatmap:
    input:
        kallisto_path= expand("results/kallisto/{unit.sample}-{unit.unit}", unit=units.itertuples()),
        sample_names="results/kallisto/"
    output:
        png="results/heatmaps/heatmap.png",
        topvar_genes="results/heatmaps/topvar_genes.tsv"
    log:
        "results/logs/heatmaps/heatmap.log",
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/get_heatmap.R"


rule get_heatmap_for_topvar_genes:
    input:
        var_genes="results/heatmaps/topvar_genes.tsv"
    output:
        var_png="results/heatmaps/top_var_genes_heatmap.png",
    log:
        "results/logs/heatmaps/top_var_genes_heatmap.log",
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/get_topvar_genes_heatmap.R"


rule get_aligned_pos:
    input:
        bam_file="results/kallisto/{sample}-{unit}/pseudoalignments.bam",
        #expand("results/kallisto/{unit.sample}-{unit.unit}/pseudoalignments.bam", unit=units.itertuples()),
    output:
        "results/QC/{sample}-{unit}.aligned.txt",
        #"results/QC/{sample}-{unit}.aligned.txt"
    log:
        "results/logs/QC/{sample}-{unit}.aligned.log",
    shell:
        "samtools view {input.bam_file} | cut -f1,3,4  > {output} 2> {log}"



rule get_read_dist:
    input:
        aligned_file="results/QC/{sample}-{unit}.aligned.txt",
    output:
        histogram="results/QC/{sample}-{unit}.histogram.html",
    log:
        "results/logs/QC/{sample}-{unit}.histogram.log",
    conda:
        "../envs/QC.yaml"
    script:
        "../scripts/get_histogram.py"
