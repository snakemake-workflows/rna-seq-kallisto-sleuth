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


use rule get_selected_transcripts_aligned_read_bins as get_aligned_read_bins with:
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


if is_3prime_experiment and config["experiment"]["3-prime-rna-seq"]["plot-qc"] != "all":

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
                "results/plots/QC/3prime-ind-QC-plot.{ind_transcripts}.html",
                category="QC",
                subcategory="global",
                caption="../report/plot-3prime-QC-histogram.rst",
                labels={
                    "plot": "3-prime read positioning",
                    "subset": "{ind_transcripts}",
                },
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
                labels={"QC-plot": "{ind_transcripts}-QC-plot"},
            ),
        params:
            each_transcript="{ind_transcripts}",
        log:
            "results/logs/QC/3prime-QC-plot.{ind_transcripts}.log",
        conda:
            "../envs/QC.yaml"
        script:
            "../scripts/plot-3prime-qc-histogram.py"
