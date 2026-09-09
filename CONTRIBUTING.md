# Contributing

For setup and checks, start with the [development guide](docs/DEVELOPMENT.md).

## Report a problem

Use the [bug report form](https://github.com/PaplomatasP/scGenesFinder/issues/new?template=bug_report.yml). Include the method and settings, expected behavior, actual behavior, console error and `sessionInfo()`. A small example that reproduces the problem is more useful than an image of the error alone.

For a proposed change, describe the analysis task it supports and how its behavior could be checked. Use the [feature request form](https://github.com/PaplomatasP/scGenesFinder/issues/new?template=feature_request.yml).

## Submit a code change

Edit the canonical application under `scgenes/` and synchronize compatibility copies with `Rscript tools/sync_app.R`. Run the applicable checks listed in the development guide. Describe the observed problem, the resulting behavior and the validation in the pull request.

Keep research outputs tied to their input and settings. If a change alters a recorded result, document the difference and update the relevant scripts, figures and checksums together. Preserve the archived function snapshot for the existing CK-p25 analysis.

The existing [AGPL v3 license](LICENSE) applies to application code. Attribute third-party code, data and annotation resources to their sources.
