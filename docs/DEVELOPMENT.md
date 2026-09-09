# Development guide

[Repository overview](../README.md) | [User guide](GETTING_STARTED.md) | [Contributing](../CONTRIBUTING.md)

## Source layout

| Path | Purpose |
| --- | --- |
| `scgenes/` | Canonical Shiny application, dependencies, assets and regression tests |
| `case-studies/ck-p25/` | Recorded biological example, input, function snapshot and results |
| `tools/` | Compatibility synchronization and check runners |
| `docs/` | User and contributor documentation |
| `images/` | Historical screenshots from an earlier interface |

Root-level `app.R`, `ui.R`, `Scripts/`, `data/` and `www/` are compatibility copies retained for existing paths. Edit the canonical files under `scgenes/`, then synchronize them from the repository root:

```sh
Rscript tools/sync_app.R
Rscript tools/check_app.R
```

The check parses the application's R files and compares the compatibility copies. Git attributes normalize application text while preserving the exact bytes of distributed case-study artifacts. The research snapshot under `case-studies/ck-p25/code_snapshot/` is intentionally frozen for the recorded analysis.

## Regression checks

Install the test dependencies in R:

```r
install.packages(c(
  "shiny", "DT", "caret", "ggplot2", "gridExtra", "ggiraph",
  "foreach", "doParallel", "zip", "png"
))
```

Then run from the repository root:

```sh
Rscript tools/check_app.R
Rscript tools/run_tests.R
Rscript case-studies/ck-p25/verify_results.R
```

The eight application scripts cover input validation, upload preview, analysis caching and errors, Run/Stop state, selected-gene k-NN classification and download settings. They use small fixtures or substitutes where needed. They do not exercise every biological method or live annotation service.

The case-study check compares saved rankings with serialized outputs, verifies plotted values and recomputes all six KEGG analyses. The separate provenance script matches the input against the public GEO matrix. Read the [case-study instructions](../case-studies/ck-p25/README.md) before regenerating saved research outputs.

## Docker

From the repository root:

```sh
docker build -t scgenes .
docker run --rm -p 8447:3838 scgenes
```

Open `http://localhost:8447` after Shiny Server starts. The image installs dependencies from `scgenes/Install_Packages/` and copies the canonical application into `/srv/shiny-server/`. The current image inherits the upstream `rocker/shiny-verse` tag and is not a frozen reproduction environment.

GitHub Actions runs the regression checks, not a complete Docker build. The legacy `update_app.sh` is a server-side deployment helper and is not executed by this workflow. Publishing a commit does not redeploy an existing application server.

## Reporting software use

Record the repository URL, commit SHA, R and package versions, input provenance and analysis settings. Cite the selection method and annotation resources used in the analysis. The repository does not assign a publication DOI or imply that the case study validates a clinical biomarker.
