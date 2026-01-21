# Replicating Violence

Replication and sensitivity analyses for Baten, Benati, and Soltysiak (2023), "Violence trends in the ancient Middle East between 12,000 and 400 BCE" (Nature Human Behaviour). This repo reproduces the paper's Model 3 (Stata intreg/Tobit), explores censoring threshold sensitivity, compares to a beta-binomial alternative, and runs robustness checks.

## Data and sources
- `S1_file_combined.xlsx` - primary analysis dataset (site-period level).
- `s41562-023-01700-y.pdf` - paper PDF.
- `0redo_viol230707nhb.do` - original Stata replication script from the authors.

## Project map
- `Model_1.R` - Model 1 replication in R (custom MLE for Stata intreg-style Tobit).
- `01_Computational_Replication.Rmd` - Model 1 replication with `survreg` and clustered SEs.
- `Model_3.R` - Model 3 replication (custom MLE, analytic weights, cluster-robust SEs) plus diagnostics and plots.
- `Model_3_Sensitivity_Threshold.R` - re-runs Model 3 across censoring points and summarizes fit and coefficient sensitivity.
- `Binomial_NEW.Rmd` - Tobit vs beta-binomial comparison with diagnostics and common-scale predictive checks.
- `Binomial.R` - older all-in-one script kept for reference.
- `Robustness check_figure2.R` - reproduces Figure 2 under alternative exclusion rules.
- `Robustness check_merged data.R` - merges duplicate rows and compares Model 3 estimates.
- `Report_Sensitivity.qmd` - main narrative report that sources Model 3 and sensitivity analyses.
- `Report_Sensitivity.html` - rendered report output.
- `Binomial_NEW.html` and `Binomial_NEW.pdf` - rendered outputs for the binomial comparison.
- `renv/` and `renv.lock` - reproducible R environment.
- `Replicating-Violence.Rproj` - RStudio project file.

## Environment setup (renv)
1) Install R 4.5.2 (see `renv.lock`).
2) Restore packages:
```bash
R -e "renv::restore()"
```
This pulls packages from the Posit CRAN mirror in `renv.lock`.

## Reproduce the analysis
Run scripts from the repo root. Most scripts print results to the console and open plots in the active graphics device.

### Model 3 replication (Stata intreg-style)
```bash
Rscript Model_3.R
```

### Censoring threshold sensitivity (Model 3)
```bash
Rscript Model_3_Sensitivity_Threshold.R
```
This script sources `Model_3.R` internally and produces fit tables and the coefficient sensitivity plot.

### Tobit vs beta-binomial comparison
```bash
R -e "rmarkdown::render('Binomial_NEW.Rmd')"
```
Rendering to PDF requires a LaTeX installation.

### Full report (Quarto)
```bash
quarto render Report_Sensitivity.qmd
```

### Robustness checks
```bash
Rscript "Robustness check_figure2.R"
Rscript "Robustness check_merged data.R"
```

### Model 1 replication (optional)
```bash
R -e "rmarkdown::render('01_Computational_Replication.Rmd')"
Rscript Model_1.R
```

## Outputs
- `Report_Sensitivity.html` - main replication and sensitivity report.
- `Binomial_NEW.html` and `Binomial_NEW.pdf` - model comparison report.

## Notes
- Scripts rely on `S1_file_combined.xlsx` in the project root.
- The Stata-intreg replication uses left-censoring at a log scale threshold (baseline `c = -3`) and analytic weights normalized to sum to N, matching the paper's implementation.
