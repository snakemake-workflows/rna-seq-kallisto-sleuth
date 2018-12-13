from snakemake.utils import validate
import pandas as pd

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


####### helpers ###########

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastqs(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_bootstrap_plots(wildcards):
    transcripts = set()
    for model in config["diffexp"]["models"]:
        results = pd.read_table(checkpoints.sleuth_diffexp.get(model=model).output[0])
        transcripts.update(results[results.qval <= config["diffexp"]["FDR"]].target_id)
    print(transcripts)
    return ["plots/bootstrap/{transcript}.bootstrap.svg".format(transcript=t)
            for t in transcripts]
