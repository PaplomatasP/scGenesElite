# Publication figures

Run an analysis, then open **Publication figures** in Run Analysis.

1. Set **Number of ranked genes** to 2–50 (default: 20). The regular ranking and classifier also support larger selections; publication panels require at most 50 for readability.
2. Choose **log2(x + 1)** for non-negative expression abundances or **As supplied** for data already transformed/scaled. The heatmap standardizes each gene across cells, then averages within groups. Constant genes are displayed as zero and listed in the export README.
3. Group by uploaded labels, or upload a CSV with `cell,sample` columns. Cell identifiers must match the expression row names, without duplicates or missing analyzed cells. Sample order follows first occurrence in the metadata CSV. Extra, unused cells are ignored. No sample identity is inferred from cell names.
4. Choose up to three ranked genes for the distribution panel. The default is the first three; change this to inspect a particular biological question.
5. Click **Generate publication figure** to preview the composite. Changing data or figure settings clears the old preview. Choose a panel and **PDF (vector)** or **PNG**, then download. The existing DPI selector controls 300/600 dpi PNG exports (default: 600).

**Download results** exports a ZIP containing the four individual plots and the composite in PDF and PNG, the full ranking, filtered data, exact plotted values, group sizes, confusion counts, predictions, classification metrics, split identifiers and a methods/session README. Plots and exports reuse the same cached classification fit. Changing only grouping, inspected genes, transformation or export settings does not train another model.

The classifier seed defaults to 20260908 and is recorded in the ZIP. Its 80% training partition, scaling and training-fold tuning follow the existing k-NN workflow. Selection occurs before the cell-level split. This is a within-dataset result, not validation on independent samples. Small classes that cannot populate a test split receive an explicit unavailable message; expression panels remain exportable.

## Reproduce Figure 2

Use the recorded CK-p25 example and SCMarker settings described in the manuscript. Set 20 genes, log2 display, seed 20260908, and upload `case-studies/ck-p25/publication_metadata.csv`. Select Il6ra, Sall1 and Lcp1. The metadata supplies descriptive sample names in the order used in the paper.

The publication export was checked against the recorded case-study heatmap values and confusion counts. The app always computes plots from the current result; it does not substitute this saved example.

## Implementation

- `Scripts/PublicationFigures.R`: validated plot data, shared plots, cached Shiny outputs and exports.
- `Scripts/PublicationUI.R`: controls and result cards.
- `Scripts/KnnClassifier.R`: numerical fit separated from rendering; split and predictions are retained.
- `Scripts/MainSelectionFun.R`: records original/mapped counts and mapped expression for both upload paths.

The old cell-level clustering/SingleR heatmap controls have been replaced by the group-mean publication heatmap. Enrichment and network tools remain in their existing sections. No new required package was added: raster exports use ragg when available and base PNG otherwise.
