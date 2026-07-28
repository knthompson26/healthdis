# Simulated Data & Code Review 

**Social connection can mitigate depression risk**

## Purpose

The original analysis in this manuscript uses restricted-access Add Health data that cannot be shared with journal reviewers. 
These scripts allow the user to run the full analysis pipeline — data simulation, bootstrap mediation analysis, and figure generation — on fabricated data that has the same structure (variable names, types, and coding) as the real analysis-ready datasets, without any real participant data.

**No real Add Health data is included or reproduced anywhere here.** 

Every value in the simulated datasets is drawn from random-number generators. 
Relationships between variables were loosely, arbitrarily calibrated only so the bootstrap mediation code has something non-null to estimate — they are **not fitted to the real data and do not represent the study's actual findings**. 
The simulated results and figures should not be interpreted as a preview of the paper's results; their only purpose is to demonstrate that the code runs correctly end-to-end.

## Contents

| File | Role |
|---|---|
| `simulate_data.R` | Generates the four fabricated analysis-ready datasets |
| `pgs.R` | Bootstrap mediation analysis: genetic risk (polygenic score) |
| `race.R` | Bootstrap mediation analysis: race |
| `ses.R` | Bootstrap mediation analysis: socio-economic position |
| `sex.R` | Bootstrap mediation analysis: sex |
| `fig2&3.rmd` | Reads the bootstrap results and produces Figures 2 and 3 |

`pgs.R`, `race.R`, `ses.R`, `sex.R`, and `fig2&3.rmd` are reviewer copies of the scripts used for the real analysis, with the minimal set of changes described in section 4 below needed to run them on simulated data.
No statistical or analytical logic was changed.

## 1. How the data was simulated

`simulate_data.R` fabricates four datasets that mirror the structure of the real analysis-ready objects produced by the study's data-preparation pipeline (`dataprep.Rmd`) and factor-score extraction (`efa.R`) that are on GitHub:

- `dat_pheno_W2`, `dat_pheno_W4` — phenotype data at Wave 2 and Wave 4 (N = 4,000 / 4,500 simulated respondents; real data N = 9,987 / 11,139)
- `dat_pheno_pgs_W2`, `dat_pheno_pgs_W4` — the above, linked to simulated genetic data for a ~45% / ~56% subset (matching the real linkage rate), including ancestry, polygenic score, genetic principal components, and an ancestry-adjusted polygenic score

Each dataset reproduces the same columns as the real data, with the same coding logic:

- Exposures: sex (0/1), race (with matching pairwise and full-sample dummy codings, including the intentional missingness pattern used for pairwise comparisons), socio-economic position (continuous and tertile dummies), and ancestry-adjusted polygenic score (continuous and tertile dummies)
- Mediators: six social connection factor scores (`family`, `school`, `friends`, `neigh`, `religion`, `social`) and their squared terms, drawn with a small, arbitrary dependence on the exposures
- Outcome: a standardized depression symptom score at each wave, generated from a synthetic linear combination of the exposures, mediators, and confounders, plus random noise
- Confounders: parental happiness, parental drinking, child health, and violence exposure

Full detail on every simulated variable, and the reasoning behind it, is documented in comments throughout `simulate_data.R`.

## 2. Setup




All five scripts point to the same local folder for reading/writing data, results, and plots:

```
/Users/qtnzknt/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Documents/healthdis/NatHB/code_review/simulated_data_scripts/
├── simulate_data.R
├── pgs.R
├── race.R
├── ses.R
├── sex.R
├── fig2&3.rmd
├── data/                              <- created by simulate_data.R
│   ├── dat_pheno_W2.RData
│   ├── dat_pheno_W4.RData
│   ├── dat_pheno_pgs_W2.RData
│   ├── dat_pheno_pgs_W4.RData
│   └── confound/all/{pgs,race,ses,sex}/   <- created by simulate_data.R, populated by the four analysis scripts
└── plots/                             <- created by fig2&3.rmd
    ├── Figure_2.png
    ├── Figure_3_PB_scaled.png
    └── Figure_3_PA.png
```

**You will need to update the directories `data_dir` to where you have saved the files.** At the moment the directories reflect my laptop paths. 

Update the `data_dir` path near the top of `simulate_data.R`, and the matching `data_dir` / `results_dir` / `plots_dir` path near the top of `pgs.R`, `race.R`, `ses.R`, `sex.R`, and `fig2&3.rmd` (each is marked `# EDITED LINE`), so that all six scripts point to the same folder.

**Required R packages:** `dplyr`, `magrittr`, `purrr`, `tidyr`, `boot`, `tidyverse`, `rlang`, `knitr`, `viridis`, `RColorBrewer`, `patchwork`, `ggplot2`, `ggpattern`, `cowplot`, `grid`, `scales`.

## 3. Running the scripts

Run in this order:

1. **`simulate_data.R`** — run once. Creates the `data/` folder and the four `.RData` files, plus the `data/confound/all/` output subfolders the analysis scripts write into.

2. **`pgs.R`, `race.R`, `ses.R`, `sex.R`** — each can be run independently, in any order. Each loads the relevant simulated dataset(s), runs the bootstrap mediation analysis across the six social-connection mediators, and saves one `.rds` results file per wave/exposure combination into `data/confound/all/<pgs|race|ses|sex>/`.

   These scripts are set up in a **quick-run mode** for review purposes: `boot_R <- 20` bootstrap resamples and `mc_reps <- 50` Monte Carlo draws per resample, instead of the `R = 1000` / `1000`-fold expansion used for the real published analysis. The full settings are extremely computationally intensive (designed to run on an HPC cluster, potentially for hours) and are unnecessary just to confirm the code runs correctly. To reproduce the full analysis settings, edit the `boot_R` and `mc_reps` values near the top of each script (see the comment directly above them).

3. **`fig2&3.rmd`** — knit after all four analysis scripts have finished. Reads every `.rds` result file from `data/confound/all/`, combines them, and produces `Figure_2.png`, `Figure_3_PB_scaled.png`, and `Figure_3_PA.png` in `plots/`.

Because of the quick-run bootstrap settings, the resulting figures will look noisier than — and may not resemble — the published results. This is expected and not a bug; it reflects the reduced number of bootstrap/Monte Carlo iterations, not an error in the code.

## 4. What was changed from the original scripts

Every line that was added or changed relative to the scripts used for the real analysis is marked in-line with the comment `# EDITED LINE` (or `# EDITED CHUNK` / `<!-- EDITED LINE -->` for whole new chunks or markdown text). Searching each file for `EDITED LINE` gives a complete, line-by-line list of every deviation. No statistical model, formula, or logic was altered — only the changes described below.

**`pgs.R`, `race.R`, `ses.R`, `sex.R`:**

- **File paths.** `load()` calls at the top and `saveRDS()` calls at the bottom were changed from Katie's private HPC cluster paths to the local `data_dir` described in section 2.
- **Bootstrap settings.** Added `boot_R` and `mc_reps` variables (see section 3) in place of the hardcoded `R = 1000` and `each = 1000` Monte Carlo expansion, so reviewers can smoke-test the code quickly.

**`fig2&3.rmd`:**

- **File paths.** Added a `results_dir` and `plots_dir` at the top, and every `readRDS()` / `ggsave()` call was updated to use them instead of Katie's private paths.
- **Axis limits.** The y-axis ranges in Figures 2 and 3 were originally hardcoded (e.g. `c(0, 0.4)`, `c(-20, 50)`) and tuned to the real published effect sizes. Because the simulated data's effect sizes are arbitrary, some results (particularly SES: Low and PGS: High) fell outside those hardcoded windows and were silently dropped from the plots. The axis limits are now computed dynamically from the data itself, so every risk factor is shown regardless of the underlying data.

**`simulate_data.R`** is a new script written for this review package; it has no "original" counterpart to diff against.

## 5. Notes and caveats

- The simulated data is entirely fabricated and contains no real participant information. Sample sizes, effect sizes, and relationships between variables are arbitrary and do not represent the study's real data or findings.
- The quick-run bootstrap settings trade statistical precision for speed. Confidence intervals and point estimates from this package will differ from — and should not be compared to — the manuscript results.
- Native American vs. White race comparisons (`natam`) are included in `race.R`'s analysis but have no case in the risk-factor label recoding in `fig2&3.rmd`'s "Join all" step, so they do not currently appear on Figures 2 or 3. This mirrors the behaviour of the original plotting script and was not something introduced by these edits; flagging it here in case it should be addressed separately.

For any questions, contact Katherine N Thompson at k.n.thompson@ucl.ac.uk. 
