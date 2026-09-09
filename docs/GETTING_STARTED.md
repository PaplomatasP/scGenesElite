# User guide

[Repository overview](../README.md) | [Case study](../case-studies/ck-p25/README.md) | [Development](DEVELOPMENT.md)

## Install and launch

Clone the repository and enter the application directory:

```sh
git clone https://github.com/PaplomatasP/scGenesFinder.git
cd scGenesFinder/scgenes
Rscript Install_Packages/install_all.R
```

The installer uses CRAN, Bioconductor and upstream research-package repositories. Packages that already load are skipped. Some dependencies require compilation, so install the toolchain appropriate to your R version, such as Rtools on Windows. Internet access is required for installation.

Start an R session in `scgenes/`, or set that as the working directory in RStudio, then run:

```r
shiny::runApp(".", host = "127.0.0.1", port = 3572, launch.browser = TRUE)
```

The browser opens the local application. The recorded CK-p25 analysis used R 4.6.1 on Windows 11; its loaded package versions are in [sessionInfo.txt](../case-studies/ck-p25/sessionInfo.txt). The regression workflow also runs on Linux. These checks do not constitute testing of every method on every R version.

## Prepare the input

Use a data frame with cells as rows, numeric gene-expression columns and class labels in the last column. A CSV can include cell identifiers in its first column. An RDS should contain the corresponding data frame, rather than a complete Seurat object.

This CSV illustrates the layout only. Its values are invented for the format example and are not study data.

```csv
cell_id,GeneA,GeneB,Labels
cell_01,2.1,0.0,control
cell_02,1.7,0.3,control
cell_03,0.5,2.8,case
cell_04,0.2,3.1,case
```

Select the matching organism, human or mouse, and gene identifier type: symbols, Ensembl gene IDs or Entrez IDs. Review genes, labels and orientation in the preview before starting selection. Use the CSV separator, header and quote controls when needed.

For a recorded biological example, download [input_ExampleData.csv](../case-studies/ck-p25/input_ExampleData.csv). This CK-p25 subset differs from the older sample under `scgenes/data/`.

## Run the CK-p25 example in the application

1. Upload the case-study CSV from the sidebar and inspect its preview.
2. Select Mouse and gene symbols. Choose SCMarker, set `geneK` and `cellK` to 10, and leave the optional variance filter and normalization off because this input already contains FPKM values.
3. Run selection and inspect the ranked genes. Use the gene-count control to change the displayed subset.
4. Export the selected data and available plots. Download bundles support 300 or 600 dpi settings.

The [scripted case study](../case-studies/ck-p25/README.md#reproduce-the-analysis) records the exact seed, function snapshot and package environment for reproduction. Its current four-panel figure uses the application's Publication figures helpers. Open **Publication figures**, select 20 genes, use log2(x + 1), upload the case-study `publication_metadata.csv`, and select Il6ra, Sall1 and Lcp1. The classifier seed is 20260908. See the [publication export guide](PUBLICATION_FIGURES.md) for downloads and interpretation.

## Interpreting results

Selection scores depend on the method. SCMarker scores are neither fold changes nor p-values for a disease comparison. Inspect expression alongside the ranking and keep the input composition with any reported gene list.

Classification uses selected genes with an 80% training split, centering and scaling, and up to five cross-validation folds within training. Gene selection precedes that split. Reported metrics describe the supplied cells and cannot establish predictive performance on independent biological samples.

Functional annotation depends on identifiers, organism, library release and background. Enrichr queries, KEGG maps and reference-based annotations can require external services or downloads. The case study's offline KEGG audit uses the mapped input genes as its background, so its values are not the application's default Enrichr results.

## Run controls and troubleshooting

Each analysis result is cached within its Shiny session. A handled selection error appears as a notification and can be retried after resolving its cause.

The Stop control is processed when the active R computation returns. It discards that run's results; it cannot interrupt a computation already occupying R.

| Symptom | Check |
| --- | --- |
| Empty preview | Confirm upload completion, file format, CSV settings and the final label column. Read any validation message. |
| Genes disappear during annotation | Check species and identifier type. Names absent from the annotation are excluded. |
| A package cannot load | Restart R. For a repeated corrupt lazy-load database error, reinstall the package named in the console. |
| Enrichment or a reference query fails | Check network access, the selected organism and the external service response. |

For a reproducible problem, [open a bug report](https://github.com/PaplomatasP/scGenesFinder/issues/new?template=bug_report.yml) with the method, settings, error text and `sessionInfo()`.
