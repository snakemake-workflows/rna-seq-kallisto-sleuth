rule kallisto_cds_index:
    input:
        "resources/transcriptome.cdna.fasta",
    output:
        "results/kallisto_cds/transcripts.idx",
    log:
        "results/logs/kallisto_cds/index.log",
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"


rule kallisto_3prime_index:
    input:
        "resources/transcriptome.3prime.fasta",
    output:
        "results/kallisto_3prime/transcripts.idx",
    log:
        "results/logs/kallisto_3prime/index.log",
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule kallisto_cds_quant:
    input:
        fq=get_trimmed,
        idx="results/kallisto_cds/transcripts.idx",
    output:
        kallisto_folder=directory("results/kallisto_cds/{sample}-{unit}"),
        bam_file="results/kallisto_cds/{sample}-{unit}/pseudoalignments.bam",
    log:
        "results/logs/kallisto_cds/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output.kallisto_folder} {params.extra} {input.fq} 2> {log}"


rule kallisto_3prime_quant:
    input:
        fq=get_trimmed,
        idx="results/kallisto_3prime/transcripts.idx",
    output:
        kallisto_folder=directory("results/kallisto_3prime/{sample}-{unit}"),
    log:
        "results/logs/kallisto_3prime/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output.kallisto_folder} {params.extra} {input.fq} 2> {log}"



rule get_heatmap:
    input:
        kallisto_path=expand("results/kallisto_3prime/{unit.sample}-{unit.unit}", unit=units.itertuples()),
        #sample_names="results/kallisto_3prime/",
    output:
        png="results/heatmaps/heatmap.png",
        matrix_file="results/heatmaps/kallisto_genematrix_file.tsv",
    params:
        sample_names="results/kallisto_3prime/",
    log:
        "results/logs/heatmaps/heatmap.log",
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/get_heatmap.R"


rule get_heatmap_for_predefine_genes:
    input:
        kallisto_genematrix_file="results/heatmaps/kallisto_genematrix_file.tsv",
        predef_genelist="resources/selected_gene_from_ref.txt",
    output:
        predefgene_png="results/heatmaps/predefgenes_heatmap.png",
    log:
        "results/logs/heatmaps/predefgenes_heatmap.log",
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/get_predefgenes_heatmap.R"


rule get_aligned_pos:
    input:
        bam_file="results/kallisto_cds/{sample}-{unit}/pseudoalignments.bam",
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
        bam_path="results/kallisto_cds/{sample}-{unit}/"
    output:
        histogram="results/QC/{sample}-{unit}.histogram.html",
    log:
        "results/logs/QC/{sample}-{unit}.histogram.log",
    conda:
        "../envs/QC.yaml"
    script:
        "../scripts/get_histogram.py"
