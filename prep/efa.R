library(polycor)
library(psych)
library(tidyverse)

# Load data
dat_fa_W2 <- readRDS("/project/rche/data/datasets/addhealth/scratch/kthompson/healthdis/efa/W2/dat_fa_W2.rds")
dat_fa_W4 <- readRDS("/project/rche/data/datasets/addhealth/scratch/kthompson/healthdis/efa/W4/dat_fa_W4.rds")

# Compute polychoric correlation matrices
polycor_matrix_W2 <- polycor::hetcor(dat_fa_W2, use = "pairwise.complete.obs", digits = 2)
saveRDS(polycor_matrix_W2, file = "/project/rche/data/datasets/addhealth/scratch/kthompson/healthdis/efa/W2/polycor_matrix_W2.rds")

polycor_matrix_W4 <- polycor::hetcor(dat_fa_W4, use = "pairwise.complete.obs", digits = 2)
saveRDS(polycor_matrix_W4, file = "/project/rche/data/datasets/addhealth/scratch/kthompson/healthdis/efa/W4/polycor_matrix_W4.rds")

# Fit and save EFA models for 1–10 factors (W2)
for (n in 1:10) {
  efa_fit_W2 <- psych::fa(
    r = polycor_matrix_W2$correlations,
    nfactors = n,
    n.obs = max(polycor_matrix_W2$n),
    n.iter = 20,
    rotate = "oblimin",
    fm = "wls"
  )
  file_path_W2 <- sprintf("/project/rche/data/datasets/addhealth/scratch/kthompson/healthdis/efa/W2/polycor_W2.efa.fit%d.rds", n)
  saveRDS(efa_fit_W2, file = file_path_W2)
  cat("Saved polychor_W2 model for nfactors =", n, "to", file_path_W2, "\n")
}

# Fit and save EFA models for 1–10 factors (W4)
for (n in 1:10) {
  efa_fit_W4 <- psych::fa(
    r = polycor_matrix_W4$correlations,
    nfactors = n,
    n.obs = max(polycor_matrix_W4$n),
    n.iter = 20,
    rotate = "oblimin",
    fm = "wls"
  )
  file_path_W4 <- sprintf("/project/rche/data/datasets/addhealth/scratch/kthompson/healthdis/efa/W4/polycor_W4.efa.fit%d.rds", n)
  saveRDS(efa_fit_W4, file = file_path_W4)
  cat("Saved polychor_W4 model for nfactors =", n, "to", file_path_W4, "\n")
}
