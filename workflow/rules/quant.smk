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


rule bwa_index:
    input:
        "resources/transcriptome_clean.cdna.fasta"
        if is_3prime_experiment
        else "resources/transcriptome.cdna.fasta",
    output:
        idx=multiext("resources/transcriptome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/transcriptome.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.17.2/bio/bwa/index"


rule bwa_mem:
    input:
        reads=get_trimmed,
        idx=multiext("resources/transcriptome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/mapped_mem/{sample}-{unit}.bam",
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.17.2/bio/bwa/mem"


rule get_mapped_canonical_transcripts:
    input:
        mapped_bam="results/mapped_mem/{sample}-{unit}.bam",
        canonical_ids="resources/canonical_ids.csv",
    output:
        canonical_mapped_bam=temp(
            "results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam"
        ),
    log:
        "results/logs/canonical_mapped_bam/{sample}-{unit}.canonical-mapped-bam.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -h -F 4 {input.mapped_bam} |  cut -f1-12 | grep -f {input.canonical_ids} | samtools view -o {output.canonical_mapped_bam}  2> {log}"


rule get_mapped_canonical_positions:
    input:
        canonical_mapped_bam="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam",
    output:
        canonical_mapped_pos=temp(
            "results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.position.txt"
        ),
    log:
        "results/logs/canonical_mapped_bam/{sample}-{unit}.canonical-mapped-pos.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {input.canonical_mapped_bam} | cut -f1,3,4,10,11  > {output}  2> {log}"


rule bwa_samtools_sort:
    input:
        "results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam",
    output:
        temp("results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.sorted.bam"),
    log:
        "results/logs/QC/{sample}-{unit}.sorted.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v1.18.3/bio/samtools/sort"


rule bwa_samtools_index:
    input:
        "results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.sorted.bam",
    output:
        temp(
            "results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.sorted.bam.bai"
        ),
    log:
        "results/logs/QC/{sample}-{unit}.sorted.index.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.18.3/bio/samtools/index"


rule get_closest_3prime_aligned_pos:
    input:
        canonical_mapped_bam="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.sorted.bam",
        canonical_mapped_bam_index="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.sorted.bam.bai",
        canonical_mapped_pos="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.position.txt",
    output:
        canonical_mapped_3prime_pos=temp(
            "results/mapped_3prime_bam/{sample}-{unit}.canonical.mapped.3prime_pos.txt"
        ),
    log:
        "results/logs/mapped_3prime_bam/{sample}-{unit}.mapped.pos.log",
    conda:
        "../envs/QC.yaml"
    script:
        "../scripts/get-3prime-max-positions.py"


rule get_closest_3prime_aligned_pos_bam:
    input:
        canonical_mapped_bam="results/canonical_mapped_bam/{sample}-{unit}.canonical.mapped.bam",
        canonical_mapped_3prime_pos="results/mapped_3prime_bam/{sample}-{unit}.canonical.mapped.3prime_pos.txt",
    output:
        canonical_mapped_3prime_bam=temp(
            "results/canonical_3prime_mapped_bam/{sample}-{unit}.canonical.3prime_mapped.bam"
        ),
    log:
        "results/logs/canonical_3prime_mapped_bam/{sample}-{unit}.canonical.3prime_mapped.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -R {input.canonical_mapped_3prime_pos} {input.canonical_mapped_bam} -o {output.canonical_mapped_3prime_bam}  2> {log}"


rule get_canonical_fastq:
    input:
        canonical_3prime_mapped_bam="results/canonical_3prime_mapped_bam/{sample}-{unit}.canonical.3prime_mapped.bam",
    output:
        canonical_fastq="results/canonical_reads/{sample}-{unit}.fastq",
    log:
        "results/logs/canonical_3prime_mapped_bam/{sample}-{unit}.canonical_3prime_mapped.fastq.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools bam2fq {input.canonical_3prime_mapped_bam} > {output.canonical_fastq}  2> {log}"


rule get_aligned_pos:
    input:
        bam_file="results/kallisto_cdna/{sample}-{unit}",
    output:
        aligned_files=temp("results/QC/{sample}-{unit}.aligned.txt"),
    log:
        "results/logs/QC/{sample}-{unit}.aligned.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {input.bam_file}/pseudoalignments.bam | cut -f1,3,4,10,11  > {output} 2> {log}"


rule kallisto_samtools_sort:
    input:
        "results/kallisto_cdna/{sample}-{unit}",
    output:
        temp("results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam"),
    log:
        "results/logs/QC/{sample}-{unit}.sorted.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input}/pseudoalignments.bam > {output} 2> {log}"


rule kallisto_samtools_index:
    input:
        "results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam",
    output:
        temp(
            "results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam.bai"
        ),
    log:
        "results/logs/QC/{sample}-{unit}.sorted.index.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.18.3/bio/samtools/index"


if is_3prime_experiment and config["experiment"]["3-prime-rna-seq"]["plot-qc"] != "all":

    rule get_selected_transcripts_aligned_read_bins:
        input:
            aligned_file="results/QC/{sample}-{unit}.aligned.txt",
            samtools_sort="results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam",
            samtools_index="results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam.bai",
            read_length="results/stats/max-read-length.json",
        output:
            fwrd_allsamp_hist_fil=temp(
                "results/QC/{sample}-{unit}.{ind_transcripts}.aligned-fwd-fil-read-bins.txt"
            ),
            rev_allsamp_hist_fil=temp(
                "results/QC/{sample}-{unit}.{ind_transcripts}.aligned-rev-fil-read-bins.txt"
            ),
        params:
            each_transcript="{ind_transcripts}",
            samples="{sample}-{unit}",
        log:
            "results/logs/QC/{sample}-{unit}.{ind_transcripts}.aligned-read-bins.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/get-sample-hist-bins.py"

    rule get_selected_transcripts_sample_QC_histogram:
        input:
            fwrd_allsamp_hist_fil=expand(
                "results/QC/{unit.sample}-{unit.unit}.{ind_transcripts}.aligned-fwd-fil-read-bins.txt",
                unit=units.itertuples(),
                ind_transcripts="{ind_transcripts}",
            ),
            rev_allsamp_hist_fil=expand(
                "results/QC/{unit.sample}-{unit.unit}.{ind_transcripts}.aligned-rev-fil-read-bins.txt",
                unit=units.itertuples(),
                ind_transcripts="{ind_transcripts}",
            ),
            read_length="results/stats/max-read-length.json",
        output:
            full_sample_QC=report(
                "results/plots/QC/3prime-QC-plot.{ind_transcripts}.html",
                category="QC",
                caption="../report/plot-3prime-QC-histogram.rst",
            ),
        params:
            each_transcript="{ind_transcripts}",
        log:
            "results/logs/QC/3prime-QC-plot.{ind_transcripts}.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/plot-3prime-qc-histogram.py"

else:

    rule get_aligned_read_bins:
        input:
            aligned_file="results/QC/{sample}-{unit}.aligned.txt",
            samtools_sort="results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam",
            samtools_index="results/kallisto-bam-sorted/{sample}-{unit}-pseudoalignments.sorted.bam.bai",
            read_length="results/stats/max-read-length.json",
        output:
            fwrd_allsamp_hist=temp(
                "results/QC/{sample}-{unit}.{ind_transcripts}.aligned-fwd-full-read-bins.txt"
            ),
            fwrd_allsamp_hist_trim=temp(
                "results/QC/{sample}-{unit}.{ind_transcripts}.aligned-fwd-trim-read-bins.txt"
            ),
            rev_allsamp_hist=temp(
                "results/QC/{sample}-{unit}.{ind_transcripts}.aligned-rev-full-read-bins.txt"
            ),
            rev_allsamp_hist_trim=temp(
                "results/QC/{sample}-{unit}.{ind_transcripts}.aligned-rev-trim-read-bins.txt"
            ),
        params:
            each_transcript="{ind_transcripts}",
            samples="{sample}-{unit}",
        log:
            "results/logs/QC/{sample}-{unit}.{ind_transcripts}.aligned-read-bins.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/get-sample-hist-bins.py"

    rule get_sample_QC_histogram:
        input:
            fwrd_allsamp_hist=expand(
                "results/QC/{unit.sample}-{unit.unit}.{ind_transcripts}.aligned-fwd-full-read-bins.txt",
                unit=units.itertuples(),
                ind_transcripts="{ind_transcripts}",
            ),
            fwrd_allsamp_hist_trim=expand(
                "results/QC/{unit.sample}-{unit.unit}.{ind_transcripts}.aligned-fwd-trim-read-bins.txt",
                unit=units.itertuples(),
                ind_transcripts="{ind_transcripts}",
            ),
            rev_allsamp_hist=expand(
                "results/QC/{unit.sample}-{unit.unit}.{ind_transcripts}.aligned-rev-full-read-bins.txt",
                unit=units.itertuples(),
                ind_transcripts="{ind_transcripts}",
            ),
            rev_allsamp_hist_trim=expand(
                "results/QC/{unit.sample}-{unit.unit}.{ind_transcripts}.aligned-rev-trim-read-bins.txt",
                unit=units.itertuples(),
                ind_transcripts="{ind_transcripts}",
            ),
            read_length="results/stats/max-read-length.json",
        output:
            full_sample_QC=report(
                "results/plots/QC/3prime-QC-plot.{ind_transcripts}.html",
                category="QC",
                caption="../report/plot-3prime-QC-histogram.rst",
            ),
        params:
            each_transcript="{ind_transcripts}",
        log:
            "results/logs/QC/3prime-QC-plot.{ind_transcripts}.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/plot-3prime-qc-histogram.py"
