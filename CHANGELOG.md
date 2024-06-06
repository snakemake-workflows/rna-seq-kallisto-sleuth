# Changelog

## [2.6.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.5.5...v2.6.0) (2024-06-05)


### Features

* Allow bam input files ([#94](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/94)) ([4a1f983](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/4a1f98320ab1b5f099941f3cd62acef7f861d631))

## [2.5.4](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.5.3...v2.5.4) (2024-01-31)


### Performance Improvements

* improve spia script and report ([#83](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/83)) ([4b3cc16](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/4b3cc16ca468ff4b05de16e906306723f6f32d09))

## [2.5.3](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.5.2...v2.5.3) (2024-01-30)


### Bug Fixes

* canonical transcript mapped read extraction ([#77](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/77)) ([52b56b0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/52b56b022729dac745724cf60266e88f55359cb3))
* the above pull request also accumulated a number of other bug fixes and updates:
  * QuantSeq data now also works with standard canonical transcripts, when MANE transcripts are not available for a species
  * some environment and wrapper updates and fixes, e.g. biomart, pysam, sleuth, datavzrd
  * some overall cleanup of the QuantSeq parts of the workflow
  * proper QuantSeq testing data, [generated with a dedicated workflow](https://github.com/dlaehnemann/create-quant-seq-testing-dataset) and [hosted on Zenodo](https://zenodo.org/doi/10.5281/zenodo.10572745), which enables a quick and useful testing for the respective parts of the workflow

## [2.5.2](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.5.1...v2.5.2) (2023-09-14)


### Bug Fixes

* simpler three prime QuantSeq cutadapt setup ([#78](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/78)) ([ecc9ab7](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/ecc9ab712b94e175a2c9e7c79b365faa98a3df44))
* update samtools.yaml to latest `1.17` and update github actions ([#75](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/75)) ([0fe7948](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/0fe79485566aad4fd9856ef62ee92ab81c6e4974))


### Performance Improvements

* bump datavzrd wrapper to 2.6.0 and general bug fixes ([#80](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/80)) ([657c465](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/657c4656d6ef45e044c0a534522e3d57d225b3e5))

## [2.5.1](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.5.0...v2.5.1) (2023-06-14)


### Bug Fixes

* Updated the get-transcript-info.R file and its dependencies ([#73](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/73)) ([e44d424](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/e44d424aede76e5443f194dd42154256e0826241))

## [2.5.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.4.3...v2.5.0) (2023-05-13)


### Features

* Support for 3-prime RNA sequencing  ([#62](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/62)) ([c06c573](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/c06c57369ea50a5e6f5f0d63a40c0ef7ae33c362))


### Bug Fixes

* fix report labels ([#72](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/72)) ([76c7fd9](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/76c7fd93a59a5f7ffae2add8b63159733a7b6c5e))

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
