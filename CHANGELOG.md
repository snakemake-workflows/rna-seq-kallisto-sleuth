# Changelog

## [2.4.3](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.4.2...v2.4.3) (2023-02-06)


### Bug Fixes

* Feature/update cutadapt wrapper ([#67](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/67)) ([29d7967](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/29d7967fac57dfdd4a8acd61d75016d8d83b5a46))

## [2.4.2](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.4.1...v2.4.2) (2022-12-02)


### Bug Fixes

* fix gene-level p-value adjustment (use Benjamini-Hochberg instead of Bonferroni-Holm) ([#64](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/64)) ([6ea1682](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/6ea1682d20b0dccd021d93359f53c3fcec9c869d))

## [2.4.1](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.4.0...v2.4.1) (2022-11-04)


### Bug Fixes

* channel order for bioconductor package download ([f57044a](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/f57044a4a1a571fcf21ce7881b32f82a3fd09265))
* correct default value for representative_transcripts and check for existence of path ([#59](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/59)) ([a85b268](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/a85b268699c74f4076d68158adfe3e3717826bbf))
* fix channel order under strict priorities ([bdbfb10](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/bdbfb103e0298b2cbee6f99802236d794d0f4797))
* fix default minimum p-value in fgsea ([#61](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/61)) ([a6a857d](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/a6a857dff26d39dc16d98c93351cf5a903f15120))
* for some rules, omit software env when caching ([#63](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/63)) ([1d2e3a9](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/1d2e3a9515701645207b57532e7045e94b3b5f3b))

## [2.4.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.3.2...v2.4.0) (2022-03-29)


### Features

* adapt to fgsea updates, configure fgsea precision by minimum achievable p-value ([dcd77ca](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/dcd77ca90ead1213acd0c293d500c18c0e579222))
