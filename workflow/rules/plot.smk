rule vg2svg:
    input:
        "{prefix}.vl.json",
    output:
        "{prefix}.svg",
    log:
        "logs/vg2svg/{prefix}.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl2svg {input} {output} 2> {log}"


rule vega_volcano_plot:
    input:
        tsv="results/tables/diffexp/{model}.transcripts.diffexp.nona.tsv",
        spec="resources/vega_volcano_plot.json",
    output:
        json="results/plots/interactive/volcano/{model}.vl.json",
    params:
        model=get_model,
        sig_level_volcano=config["diffexp"]["sig-level"]["volcano-plot"],
        primary_variable=lambda wc: config["diffexp"]["models"][wc.model][
            "primary_variable"
        ],
    log:
        "logs/vega-plots/volcano/{model}.log",
    conda:
        "../envs/vega.yaml"
    script:
        "../scripts/vega_plot_volcano.py"


rule dropna:
    input:
        tsv="results/tables/diffexp/{model}.transcripts.diffexp.tsv",
    output:
        tsv="results/tables/diffexp/{model}.transcripts.diffexp.nona.tsv",
    log:
        "logs/dropna/{model}.log",
    run:
        df = pd.read_csv(input.tsv, sep="\t").dropna()
        df.to_csv(output.tsv, sep="\t", index=False)
