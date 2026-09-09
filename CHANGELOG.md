# Changes

## 2026-09-09

Added the Publication figures workspace: selection counts, group-mean expression heatmaps, individual gene distributions and cached k-NN confusion matrices. Users can supply cell-to-sample metadata, choose the expression transform and export vector PDF or 300/600 dpi PNG plots. Download bundles include plotted values, rankings, predictions, split identifiers, settings and methods notes.

The classifier records its seed and reuses the same fit for the displayed plots and downloads. CSV and RDS selection paths retain original and mapped gene counts. Publication regression checks bring the suite to eight scripts. Both supported application locations contain the same current files.

Updated the README with the current four-panel CK-p25 figure, its reproducible export script and data bundle. Updated GitHub links to the current scGenesFinder repository name.

## 2026-09-06

Updated the application with sidebar CSV/RDS upload, a server-driven preview and input validation. Selection results are cached per session, and handled analysis errors can be retried without closing the session. Run/Stop controls report readiness and prevent duplicate clicks; Stop is processed after the active R computation returns.

The k-nearest neighbours classifier now uses only the selected gene columns. It uses an 80% training split, centering and scaling, and up to five cross-validation folds within training. Confusion-matrix colors distinguish correct from incorrect classifications. Download bundles support plot resolution settings of 300 or 600 dpi.

Added regression checks for input, preview state, analysis errors, run controls, selected-gene classification and exports. Added CI and a check that root compatibility files match the canonical application under `scgenes/`. Docker uses the canonical dependency installers. Existing screenshots in `images/` are historical.

Added the CK-p25 case study with its exact input subset, original function snapshot, package versions, reproducible scripts, numerical results and vector/600 dpi figures. The public package excludes the manuscript draft and local session histories.
