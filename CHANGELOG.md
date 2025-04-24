# Changelog

## [2.11.1](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.11.0...v2.11.1) (2025-04-24)


### Bug Fixes

* Remove model filtering of logcount matrix ([#144](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/144)) ([de347fa](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/de347fa4502ce7b15c62bd3a4fa0f4525f2f86f8))

## [2.11.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.10.0...v2.11.0) (2025-04-23)


### Features

* Add TPM matrix generation and integration into workflow ([#146](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/146)) ([3c6bdac](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/3c6bdac2396b268e102a9f7dd277b4559863e381))


### Bug Fixes

* make postprocess sorting work with more than one signed_pi value column ([#150](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/150)) ([4b3417a](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/4b3417abb8aba5d99df29899d6e500c4eebe512f))
* use `verify_integrity=True` in pandas `set_index()` for samples and units to check for duplicates ([#149](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/149)) ([9657c75](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/9657c75d559bc41f5fe649ad631c59a50733bc55))

## [2.10.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.9.0...v2.10.0) (2025-03-19)


### Features

* Enrichment and pathway scatter-plots ([#136](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/136)) ([29684af](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/29684af31dfe7fb6fb05ed1b0fa52a2f356c8091))


### Bug Fixes

* move summing of go-term fold changes out of log-space ([#141](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/141)) ([66bed4a](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/66bed4af4fc061a4ec6ff5aacfea0dd1a826bfa1))

## [2.9.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.8.4...v2.9.0) (2025-03-11)


### Features

* Make PCA plots interactive and switch to HTML format ([#137](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/137)) ([d0e6e78](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/d0e6e785d0b2518f11efcf3d18a44a53bb578035))


### Bug Fixes

* add support for drosophila melanogaster and other non-standard species ([#142](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/142)) ([dedfdf7](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/dedfdf714c28c22d58546eb369536957768e378c))
* Improve error handling for missing or empty column values in dataframe with duplicate definitions ([#138](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/138)) ([7f3d590](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/7f3d59000897b5831ed56cb064789d21b95b8935))

## [2.8.4](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.8.3...v2.8.4) (2024-12-18)


### Bug Fixes

* handle NA values correctly in meta compare pathways ([#134](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/134)) ([71778df](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/71778dfacc9f662ba64916fb9a99b6107912f05f))

## [2.8.3](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.8.2...v2.8.3) (2024-12-11)


### Bug Fixes

* Update to latest datavzrd wrapper ([#126](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/126)) ([8b57d60](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/8b57d60b713fe3994f95afc069d59cb089371524))
* use latest stable datavzrd, fixing a bug with numpy values in configs; increase default ensembl release version to 113 (112 currently has biomart issues) ([#130](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/130)) ([b11f100](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/b11f10061cbc38cf693ff6d919eba4e764630ac9))

## [2.8.2](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.8.1...v2.8.2) (2024-10-18)


### Bug Fixes

* Convert pval and qval to heatmap ([#123](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/123)) ([478c759](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/478c75947e87aac8feb0d71d6b12e112bae58570))
* Move NA values to end of table ([#124](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/124)) ([df807a5](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/df807a51b312fe14b6cd2e6f9f9b0d048843c859))
* Update max in memory rows for diffexp report ([#122](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/122)) ([830e678](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/830e678fbf5a055c6d23cdb64b6b31953d7a4031))


### Performance Improvements

* Update to latest datavzrd wrapper ([#125](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/125)) ([3a4c020](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/3a4c020f9f1927cdd62c34543f48c97db7869df3))
* Update to latest datavzrd wrapper version ([#120](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/120)) ([7b2767c](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/7b2767cc5ebcbd34aea096047712a179b60b868e))

## [2.8.1](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.8.0...v2.8.1) (2024-09-20)


### Performance Improvements

* Update datavzrd wrapper version ([#118](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/118)) ([e498079](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/e49807901c2d934acef794bf93b3b277c7bbe868))

## [2.8.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.7.2...v2.8.0) (2024-09-20)


### Features

* sort datavzrd tables ([#115](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/115)) ([1ccc3fa](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/1ccc3fa626e725fa1cbc19464c0c441f74e5edb3))


### Bug Fixes

* Move signed_pi_value columns to the end (again) ([#117](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/117)) ([ac4a73b](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/ac4a73b5a4102f96aeb18edc7df523446ef35533))

## [2.7.2](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.7.1...v2.7.2) (2024-09-11)


### Bug Fixes

* Fix meta comparison and wildcard issue ([#113](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/113)) ([3f8f126](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/3f8f1265ea035038dba2e1115f46e0ad01717079))
* Fix meta comparison model and label path ([#111](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/111)) ([27b9cb1](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/27b9cb136eb5fa9ecf4d0843c57be791f56ab730))

## [2.7.1](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.7.0...v2.7.1) (2024-08-22)


### Bug Fixes

* Move signed_pi_value_* columns into detail mode that dont include the primary variable ([#110](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/110)) ([a159b4c](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/a159b4cf239488bd284ed67f383f11f36e30f056))
* Move signed_pi_value_* columns to the end and move non primary variable columns to detail mode ([#106](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/106)) ([8b2e3fe](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/8b2e3fe78a560783babb850697642de41d134b1e))
* Remove IHW outputs from report ([#107](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/107)) ([aa891fa](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/aa891faf4d372b2821fd67b62908bb516bff3eba))
* Update datavzrd wrapper version ([#109](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/109)) ([1fe8c2c](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/1fe8c2cf6734e179ce4f0da5501eb55e4de53128))

## [2.7.0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/compare/v2.6.0...v2.7.0) (2024-08-15)


### Features

* Improve datavzrd tables ([#93](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/93)) ([93512b8](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/93512b8c7ae8b2fbc1bba608168dd0309ee5e0b1))


### Bug Fixes

* Fix missing output in spia.R when no significant genes are found ([#103](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/103)) ([bc0d017](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/bc0d017a1cc67118142711a2fc4fb0bb31218fe2))
* Handle missing bam columns in units.tsv ([#105](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/105)) ([bae88d0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/bae88d06b2bb2ef606175fb231acebe9491d05cc))
* Remove non-existent outputs in spia rule ([#102](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/102)) ([0fbb930](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/0fbb93065ef16f593dfbd0eb6332eb18b9237e60))
* update to latest datavzrd ([417ec3b](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/417ec3b26a69d1549f7dfbc9b1c9f0b4d99209b7))


### Performance Improvements

* datavzrd wrapper `v3.12.1`, offer-excel configurable, free disk space for CI, dynamic sleuth_init mem_mb, pure download rules as localrules ([#92](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/92)) ([70850fb](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/70850fb7f573e1868dc9400a0af8d8ffe86435e6))
* Update datavzrd wrapper ([#98](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/98)) ([e5eb0e0](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/e5eb0e041a220901c7a7fcba60d8a963749319b9))
* Update samtools fast separate wrapper ([#100](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/issues/100)) ([65d8f41](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/commit/65d8f4132c3606bf620b8bdc1ffd2785d6f7c17e))

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
