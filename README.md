# scGenesFinder

scGenesFinder is an R Shiny application for gene selection and interpretation of single-cell RNA sequencing data. This repository retains its original GitHub name, **scGenesElite**.

The application combines differential expression methods, methods developed for single-cell data, classifier variable importance, SHAP-based selection and ensemble ranking. Selected genes can be inspected through expression plots, k-nearest neighbours classification, Enrichr queries, KEGG maps and interaction graphs.

## Run locally

The canonical application is in `scgenes/`. From the repository root:

```sh
cd scgenes
Rscript Install_Packages/install_all.R
Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3572, launch.browser=TRUE)'
```

In RStudio, set the working directory to `scgenes/` and run `shiny::runApp()` after installing dependencies. GitHub packages may require a compiler toolchain, such as Rtools on Windows. Installation uses CRAN, Bioconductor and the upstream research-package repositories listed in `Install_Packages/`. The case study was run with R 4.6.1 on Windows 11; its exact loaded package versions are recorded separately.

The root-level `app.R`, `ui.R`, `Scripts/`, `data/` and `www/` are synchronized compatibility copies for existing paths. Make application changes under `scgenes/` and run `Rscript tools/sync_app.R` before committing.

## Input and analysis

Upload CSV or RDS data from the sidebar. Use cells as rows, numeric gene-expression columns, and cell labels in the final column. A leading CSV cell-identifier column is recognized as row names. Select the correct organism, human or mouse, and identifier type. The preview displays the uploaded data before analysis.

Choose a selection method and preprocessing settings, then run the analysis. The returned ranking controls the selected gene set used by downstream outputs. Download bundles support 300 or 600 dpi plots. Each run is cached within its Shiny session, and handled selection errors are displayed without closing the session.

The Stop control discards the completed run after R can process the request. It does not interrupt an R computation already in progress. Classification splits cells after gene selection and therefore does not estimate performance on independent biological samples. Enrichment and reference-based annotation require network access.

## CK-p25 microglia case study

The [case-study directory](case-studies/ck-p25/README.md) contains the exact input subset, archived analysis functions, scripts, gene rankings, pathway-test results and publication-size figures.

| Quantity | Recorded result |
| --- | --- |
| Source | GSE103334, CK-p25 mouse hippocampal microglia |
| Input | 384 cells, 2,000 gene columns, four source samples |
| Annotation | 1,860 mapped genes |
| SCMarker output | 302 genes, an 83.8% reduction from the mapped input |
| Week-2 subset | 192 cells, 374 selected genes |
| Top-50 overlap between runs | 7 genes |

![SCMarker ranking and expression across the four CK-p25 source samples](case-studies/ck-p25/figures/scGenesFinder_case_study_600dpi.png)

This figure was prepared from the saved analysis results. It is not a browser screenshot. The heatmap shows gene-wise standardized log2(FPKM + 1) values averaged within each source sample. The separate KEGG audit found no terms with FDR below 0.05. The example documents selection and expression inspection, without establishing independent disease biomarkers or comparative predictive performance.

## Checks

From the repository root, run:

```sh
Rscript tools/check_app.R
Rscript tools/run_tests.R
Rscript case-studies/ck-p25/verify_results.R
```

The application regression checks require `shiny`, `caret`, `ggplot2`, `gridExtra`, `foreach`, `doParallel`, `zip`, `png` and their dependencies. They cover CSV/RDS validation, preview state, cached analysis errors, Run/Stop state, selected-gene classification and download settings. They do not exercise every biological method or external service. GitHub Actions runs these checks on pushes and pull requests.

## Docker

```sh
docker build -t scgenes .
docker run --rm -p 8447:3838 scgenes
```

The image installs the canonical dependencies and copies `scgenes/` into Shiny Server. A full Docker build is separate from the local regression checks. Publishing this repository does not update a running server.

## License and attribution

Application code is distributed under the existing [GNU AGPL v3 license](LICENSE). Public study data and annotation resources retain their source attribution and applicable terms. See the [case-study provenance](case-studies/ck-p25/README.md#data-provenance). Historical interface images remain under `images/`; they describe an earlier interface.