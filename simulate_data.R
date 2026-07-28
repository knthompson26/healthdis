# =============================================================================
# simulate_data.R
#
# Purpose: Generate a fully synthetic ("dummy") version of the four analytic
# datasets produced by dataprep.Rmd (with mediator factor scores from efa.R),
# so that the user can run pgs.R, race.R, ses.R, and sex.R
# without access to the restricted Add Health data.
#
# No real participant data is used anywhere in this script. Every value is
# drawn from random-number generators. Associations between variables are
# loosely calibrated (not fit to real data) purely so the downstream
# mediation/bootstrap code produces sensible, non-degenerate output -- they
# do NOT represent the study's actual findings and should not be
# interpreted as such.
#
# Output: four .RData files, each restoring one object under its original
# name when loaded:
#   dat_pheno_W2.RData      -> dat_pheno_W2
#   dat_pheno_W4.RData      -> dat_pheno_W4
#   dat_pheno_pgs_W2.RData  -> dat_pheno_pgs_W2
#   dat_pheno_pgs_W4.RData  -> dat_pheno_pgs_W4
#
# These are written to `data_dir` below, which matches the path already set
# at the top of the reviewer copies of pgs.R, race.R, ses.R, and sex.R.
# If you move the project folder, update `data_dir` here AND the `data_dir`
# line (marked "# EDITED LINE") near the top of each of those four scripts.
# =============================================================================

library(dplyr)
library(magrittr)  # for %>%, used inside add_pgs()

set.seed(20260728)  # reproducible dummy data

# ---- 0. Where to write the simulated data ----------------------------------
data_dir <- "/Users/qtnzknt/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Documents/healthdis/NatHB/code_review/simulated_data_scripts/data/"
out_dir  <- paste0(data_dir, "confound/all")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
for (sub in c("pgs", "race", "ses", "sex")) {
  dir.create(file.path(out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

# ---- 1. Sample sizes --------------------------------------------------------
# Real N (for reference, not reproduced exactly): W2 phenotype N=9,987,
# W4 phenotype N=11,139; W2 PGS-linked subset N=4,478 (~45%),
# W4 PGS-linked subset N=6,242 (~56%). We use smaller but still ample N so
# the bootstrap in each script runs quickly.
n_w2 <- 4000
n_w4 <- 4500
pgs_prop_w2 <- 0.45
pgs_prop_w4 <- 0.56

# ---- 2. Helper: simulate one wave of phenotype data -------------------------
# Reproduces (with fabricated values) every column documented in
# dataprep.Rmd's final dat_pheno_W2 / dat_pheno_W4 objects:
#   AID, H1GI1Y, SCID, RACE1, SES, CESD_W{2,4},
#   family, school, friends, neigh, religion, social (+ squared terms),
#   sex, race, UM, UY,
#   black, natam, asian, SES_cat, SESlow, SESmid,
#   black_cov, natam_cov, asian_cov,
#   happyP, drinkP_bin, health_bin, violence
make_pheno <- function(n, wave = c("W2", "W4")) {
  wave <- match.arg(wave)

  AID    <- sprintf("AID%06d", sample(100000:999999, n))
  SCID   <- sample(1:132, n, replace = TRUE)
  H1GI1Y <- sample(74:83, n, replace = TRUE)          # birth year 1974-1983 

  # sex: 0/1 numeric as required by the models (NOT a factor).
  # Coding direction (0=male/1=female or vice versa) is not documented in the
  # original scripts; adjust the label below if Katie's real coding differs.
  sex <- rbinom(n, 1, 0.51)                            # 0 = male, 1 = female

  # race: RACE1 1=White, 2=Black, 3=Native American, 4=Asian
  RACE1 <- sample(1:4, n, replace = TRUE, prob = c(0.60, 0.25, 0.03, 0.12))

  # SES: continuous composite (~SESPC_AL), bounded 0-100 for plausibility
  SES   <- pmin(pmax(rnorm(n, mean = 50, sd = 15), 0), 100)
  SES_tertile <- dplyr::ntile(SES, 3)                   # 1=Low, 2=Middle, 3=High
  SES_cat <- c("Low", "Middle", "High")[SES_tertile]

  # Wave-1 confounders
  happyP     <- rbinom(n, 1, 0.75)
  drinkP_bin <- rbinom(n, 1, 0.25)
  health_bin <- rbinom(n, 1, 0.15)
  violence   <- as.numeric(scale(rpois(n, lambda = 1)))

  # race/SES dummy coding, mirroring dataprep.Rmd's derivation logic exactly:
  # - black/natam/asian: pairwise vs. White, NA for the other minority groups
  #   (used as the focal 0/1 exposure, one at a time, in race.R)
  black <- dplyr::case_when(RACE1 == 2 ~ 1, RACE1 == 1 ~ 0, TRUE ~ NA_real_)
  natam <- dplyr::case_when(RACE1 == 3 ~ 1, RACE1 == 1 ~ 0, TRUE ~ NA_real_)
  asian <- dplyr::case_when(RACE1 == 4 ~ 1, RACE1 == 1 ~ 0, TRUE ~ NA_real_)

  # - black_cov/natam_cov/asian_cov: full-sample 0/1 dummies (no NA), used as
  #   covariates in ses.R / sex.R / pgs.R
  black_cov <- ifelse(RACE1 == 2, 1, 0)
  natam_cov <- ifelse(RACE1 == 3, 1, 0)
  asian_cov <- ifelse(RACE1 == 4, 1, 0)

  race <- factor(RACE1, levels = c(1, 2, 3, 4))

  # SESlow/SESmid: pairwise vs. High, NA for the excluded tertile
  # (used as the focal 0/1 exposure, one at a time, in ses.R)
  SESlow <- dplyr::case_when(SES_cat == "Low"    ~ 1, SES_cat == "High" ~ 0, TRUE ~ NA_real_)
  SESmid <- dplyr::case_when(SES_cat == "Middle" ~ 1, SES_cat == "High" ~ 0, TRUE ~ NA_real_)

  # mediators: social-connection latent factor scores (mimicking the lavaan
  # CFA output in dataprep.Rmd). Given a small, plausible dependence on
  # exposures so the mediation models have something non-null to estimate.
  race_penalty <- dplyr::case_when(RACE1 == 2 ~ -0.25,
                                    RACE1 == 3 ~ -0.20,
                                    RACE1 == 4 ~ -0.10,
                                    TRUE ~ 0)
  ses_boost <- 0.01 * (SES - 50)

  family   <- rnorm(n, 0.15 * sex + race_penalty + ses_boost, 1)
  school   <- rnorm(n, 0.10 * sex + race_penalty + ses_boost, 1)
  friends  <- rnorm(n, 0.05 * sex + race_penalty + ses_boost, 1)
  neigh    <- rnorm(n, race_penalty + ses_boost, 1)
  religion <- rnorm(n, ses_boost, 1)
  social   <- as.numeric(scale(family + school + friends + neigh + religion + rnorm(n, 0, 0.5)))

  family2   <- family^2
  school2   <- school^2
  friends2  <- friends^2
  neigh2    <- neigh^2
  religion2 <- religion^2
  social2   <- social^2

  UM <- rnorm(n)
  UY <- rnorm(n)

  # outcome: standardized CES-D depression score
  lin_pred <- 0.35 * sex - 0.015 * (SES - 50) +
    0.30 * black_cov + 0.15 * asian_cov + 0.10 * natam_cov -
    0.35 * social - 0.10 * family - 0.05 * school +
    0.25 * violence - 0.15 * happyP + 0.20 * drinkP_bin + 0.25 * health_bin +
    rnorm(n, 0, 1)
  CESD <- as.numeric(scale(lin_pred))

  out <- data.frame(
    AID, H1GI1Y, SCID, RACE1, SES,
    family, school, friends, neigh, religion, social,
    sex,
    family2, school2, friends2, neigh2, religion2, social2,
    race, UM, UY,
    black, natam, asian,
    SES_cat, SESlow, SESmid,
    black_cov, natam_cov, asian_cov,
    happyP, drinkP_bin, health_bin, violence,
    stringsAsFactors = FALSE
  )
  out[[paste0("CESD_", wave)]] <- CESD
  out
}

# ---- 3. Helper: link a subset of a phenotype df to fabricated PGS data ------
# Reproduces the columns added to dat_pheno_pgs_W2 / dat_pheno_pgs_W4:
#   ancestry, pgs, pgs_cat, pgs2, pgs3, PC1-PC10,
#   pgs_reg, pgs_reg_cat, pgs_reg2, pgs_reg3, pgs_reg_cat_factor
add_pgs <- function(df, prop) {
  n_total <- nrow(df)
  keep <- sort(sample.int(n_total, size = round(prop * n_total)))
  d <- df[keep, ]
  n <- nrow(d)
  rownames(d) <- NULL

  ancestry <- sample(c("EUR", "AFR", "AMR", "EAS"), n, replace = TRUE,
                      prob = c(0.55, 0.25, 0.13, 0.07))

  PCs <- as.data.frame(matrix(rnorm(n * 10, 0, 0.05), nrow = n))
  names(PCs) <- paste0("PC", 1:10)

  pgs      <- rnorm(n, 0, 1)                      # raw polygenic score
  pgs_cat  <- dplyr::ntile(pgs, 3)                 # pooled tertile (1=low,2=mid,3=high)
  pgs2 <- dplyr::case_when(pgs_cat == 2 ~ 1, pgs_cat == 1 ~ 0, TRUE ~ NA_real_)
  pgs3 <- dplyr::case_when(pgs_cat == 3 ~ 1, pgs_cat == 1 ~ 0, TRUE ~ NA_real_)

  d <- cbind(d, ancestry, pgs, pgs_cat, pgs2, pgs3, PCs)

  # ancestry-adjusted PGS: residualize pgs on the 10 PCs within each
  # ancestry group, then standardize (mirrors the real pipeline)
  d$pgs_reg <- NA_real_
  for (anc in unique(d$ancestry)) {
    idx <- d$ancestry == anc
    fit <- lm(pgs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
              data = d[idx, ])
    d$pgs_reg[idx] <- as.numeric(scale(resid(fit)))
  }

  d <- d %>%
    dplyr::group_by(ancestry) %>%
    dplyr::mutate(pgs_reg_cat = dplyr::ntile(pgs_reg, 3)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      pgs_reg2 = dplyr::case_when(pgs_reg_cat == 2 ~ 1, pgs_reg_cat == 1 ~ 0, TRUE ~ NA_real_),
      pgs_reg3 = dplyr::case_when(pgs_reg_cat == 3 ~ 1, pgs_reg_cat == 1 ~ 0, TRUE ~ NA_real_),
      pgs_reg_cat_factor = factor(pgs_reg_cat, levels = 1:3, labels = c("Low", "Middle", "High"))
    ) %>%
    as.data.frame()

  d
}

# ---- 4. Build the four datasets ---------------------------------------------
dat_pheno_W2 <- make_pheno(n_w2, "W2")
dat_pheno_W4 <- make_pheno(n_w4, "W4")

dat_pheno_pgs_W2 <- add_pgs(dat_pheno_W2, pgs_prop_w2)
dat_pheno_pgs_W4 <- add_pgs(dat_pheno_W4, pgs_prop_w4)

# add a small, plausible PGS -> depression association so pgs.R's
# mediation estimates aren't identically null
dat_pheno_pgs_W2$CESD_W2 <- as.numeric(scale(dat_pheno_pgs_W2$CESD_W2 + 0.25 * dat_pheno_pgs_W2$pgs_reg))
dat_pheno_pgs_W4$CESD_W4 <- as.numeric(scale(dat_pheno_pgs_W4$CESD_W4 + 0.25 * dat_pheno_pgs_W4$pgs_reg))

# ---- 5. Save --------------------------------------------------------------
save(dat_pheno_W2,     file = paste0(data_dir, "dat_pheno_W2.RData"))
save(dat_pheno_W4,     file = paste0(data_dir, "dat_pheno_W4.RData"))
save(dat_pheno_pgs_W2, file = paste0(data_dir, "dat_pheno_pgs_W2.RData"))
save(dat_pheno_pgs_W4, file = paste0(data_dir, "dat_pheno_pgs_W4.RData"))

cat("Simulated data written to:", data_dir, "\n")
cat("  dat_pheno_W2:    ", nrow(dat_pheno_W2), "rows x", ncol(dat_pheno_W2), "cols\n")
cat("  dat_pheno_W4:    ", nrow(dat_pheno_W4), "rows x", ncol(dat_pheno_W4), "cols\n")
cat("  dat_pheno_pgs_W2:", nrow(dat_pheno_pgs_W2), "rows x", ncol(dat_pheno_pgs_W2), "cols\n")
cat("  dat_pheno_pgs_W4:", nrow(dat_pheno_pgs_W4), "rows x", ncol(dat_pheno_pgs_W4), "cols\n")
