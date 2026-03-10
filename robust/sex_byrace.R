library(purrr)
library(tidyr)
library(boot)
library(tidyverse)
library(dplyr)
library(rlang)
library(magrittr)

set.seed(123)
remove(list = ls())

# load data
load("/home/thom1336/healthdis/data/dat_pheno_W2.RData")
load("/home/thom1336/healthdis/data/dat_pheno_W4.RData")

# list of mediators
mediators <- c("family", "school", "friends", "neigh", "religion", "social")

# create racial minority specific datasets based on dummy vars
for(i in c("black", "asian", "natam")){
  assign(paste0("dat_pheno_W2_", i), dat_pheno_W2 %>% dplyr::filter(.data[[i]] == 1))
  assign(paste0("dat_pheno_W4_", i), dat_pheno_W4 %>% dplyr::filter(.data[[i]] == 1))
}

# create white dataset
dat_pheno_W2_white <- dat_pheno_W2 %>% dplyr::filter(black == 0 & asian == 0 & natam == 0)
dat_pheno_W4_white <- dat_pheno_W4 %>% dplyr::filter(black == 0 & asian == 0 & natam == 0)

######----------------------- Runs a loop for each race-specific dataset -----------------------#######

for(race in c("black", "asian", "natam", "white")){

# list of waves and their corresponding data and outcome names
waves <- list(
  W2 = list(data = get(paste0("dat_pheno_W2_", race)), outcome = "CESD_W2"),
  W4 = list(data = get(paste0("dat_pheno_W4_", race)), outcome = "CESD_W4")
)

# function to run counterfactual calculations and Monte Carlo simulations within bootstrapping 
bootstrap_stat <- function(data, indices, wave_name, mediators) {
  
  #--------------------------------------------
  # Bootstrap sample (original sample)
  #--------------------------------------------
  dat <- data[indices, ]
  dat$original <- 1
  dat$UM <- rnorm(nrow(dat))  # mediator error (original)
  dat$UY <- rnorm(nrow(dat))  # outcome error (original)
  
  outcome_var <- waves[[wave_name]]$outcome
  
  #--------------------------------------------
  # Monte Carlo expansion (inner loop)
  #--------------------------------------------
  dat_expanded <- dat[rep(1:nrow(dat), each = 1000), ]  
  dat_expanded$original <- rep(c(1, rep(0, 999)), times = nrow(dat))
  dat_expanded$UM <- rnorm(nrow(dat_expanded))
  dat_expanded$UY <- rnorm(nrow(dat_expanded))
  
  #--------------------------------------------
  # Counterfactual loop
  #--------------------------------------------
  
  # create a list to store model results
  model_list_regsex <- list()
  
  #----------------------------------------------------
  # Step 1: Variable for social (M) when sex (X) = 0
  #----------------------------------------------------
  
  for (m in mediators) {
    # fit model
    modelM <- lm(as.formula(paste(m, "~ sex")), data = dat)  # fit on original bootstrapped sample
    
    # store model in list
    model_list_regsex[[paste0("model_", m, "_sex")]] <- modelM
    
    # extract intercept and residual SD
    intercept <- coef(modelM)[["(Intercept)"]]
    resid_sd <- summary(modelM)$sigma
    
    # simulate mediator under sex = 0 (on expanded)
    new_var <- paste0(m, "_sex0")
    dat_expanded[[new_var]] <- intercept + resid_sd * dat_expanded$UM # applied to all expanded simulations
    
    print(summary(dat[[m]]))
    print(summary(dat_expanded[[new_var]]))
  }
  
  #------------------------------------
  # Step 2: Model for Dep (Y)
  #------------------------------------
  
  walk(mediators, function(m) {
    # construct formula for the model
    formula_str <- paste0(
      outcome_var, " ~ ",
      "sex + ",                                                    # x
      m, " + ",                                                    # m
      m, ":sex + ",                                                # m*x
      m, "2 + ",                                                   # m2
      "H1GI1Y + ",                                                 # c1
      "SES + ",                                                    # c2
      "happyP + drinkP_bin + health_bin + violence + ",           
      
      # sex × covariates
      "sex:H1GI1Y + ",                                             # c1*x
      "sex:SES + ",                                                # c2*x 
      "sex:happyP + ",
      "sex:drinkP_bin + ",
      "sex:health_bin + ",
      "sex:violence + ",
      
      # coavriate x covariate
      "H1GI1Y:SES + ",                                             
      "H1GI1Y:happyP + ",
      "H1GI1Y:drinkP_bin + ",
      "H1GI1Y:health_bin + ",
      "H1GI1Y:violence + ",
      "SES:happyP + ",
      "SES:drinkP_bin + ",
      "SES:health_bin + ",
      "SES:violence + ",
      "happyP:drinkP_bin + ",
      "happyP:health_bin + ",
      "happyP:violence + ",
      "drinkP_bin:health_bin + ",
      "drinkP_bin:violence + ",
      "health_bin:violence"
    )
    
    formula_obj <- as.formula(formula_str)
    modelY <- lm(formula_obj, data = dat)  # fit outcome model on original sample
    model_nameY <- paste0("model_", outcome_var, "_sex_", m)
    assign(model_nameY, modelY, envir = .GlobalEnv)
  })
  
  #-----------------------------------------------------------------------
  # Step 2A: Variable for Y, when X=0 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_sex_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    sex0_var <- sym(paste0(m, "_sex0"))
    output_var_00 <- sym(paste0(outcome_var, "_00_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_00 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # X = sex = 0
          (coefs[["sex"]] %||% 0) * 0 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                      * !!sex0_var +
          (coefs[[paste0("sex:", m)]] %||% 0)      * 0 * !!sex0_var +
          (coefs[[paste0(m, "2")]] %||% 0)         * (!!sex0_var)^2 +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]    %||% 0) * H1GI1Y +
          (coefs[["SES"]]       %||% 0) * SES +
          (coefs[["happyP"]]    %||% 0) * happyP +
          (coefs[["drinkP_bin"]]%||% 0) * drinkP_bin +
          (coefs[["health_bin"]]%||% 0) * health_bin +
          (coefs[["violence"]]  %||% 0) * violence +
          
          # interactions with sex (X = 0 → all vanish)
          (coefs[["sex:H1GI1Y"]]    %||% 0) * H1GI1Y   * 0 +
          (coefs[["sex:SES"]]       %||% 0) * SES       * 0 +
          (coefs[["sex:happyP"]]    %||% 0) * happyP    * 0 +
          (coefs[["sex:drinkP_bin"]]%||% 0) * drinkP_bin * 0 +
          (coefs[["sex:health_bin"]]%||% 0) * health_bin * 0 +
          (coefs[["sex:violence"]]  %||% 0) * violence  * 0 +
          
          # covariate–covariate interactions
          (coefs[["H1GI1Y:SES"]]         %||% 0) * H1GI1Y * SES +
          (coefs[["H1GI1Y:happyP"]]      %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]  %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]  %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]    %||% 0) * H1GI1Y * violence +
          (coefs[["SES:happyP"]]         %||% 0) * SES * happyP +
          (coefs[["SES:drinkP_bin"]]     %||% 0) * SES * drinkP_bin +
          (coefs[["SES:health_bin"]]     %||% 0) * SES * health_bin +
          (coefs[["SES:violence"]]       %||% 0) * SES * violence +
          (coefs[["happyP:drinkP_bin"]]  %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]  %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]    %||% 0) * happyP * violence +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
          # residual
          rmse * UY
      )
  })
  
  #----------------------------------------------------------------------
  # Step 2B: Variable for Y, when X=1 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_sex_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    
    sex0_var   <- sym(paste0(m, "_sex0"))
    output_var_10 <- sym(paste0(outcome_var, "_10_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_10 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # X = sex = 1
          (coefs[["sex"]] %||% 0) * 1 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                 * !!sex0_var +
          (coefs[[paste0("sex:", m)]] %||% 0) * 1 * !!sex0_var +
          (coefs[[paste0(m, "2")]] %||% 0)    * (!!sex0_var)^2 +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]    %||% 0) * H1GI1Y +
          (coefs[["SES"]]       %||% 0) * SES +
          (coefs[["happyP"]]    %||% 0) * happyP +
          (coefs[["drinkP_bin"]]%||% 0) * drinkP_bin +
          (coefs[["health_bin"]]%||% 0) * health_bin +
          (coefs[["violence"]]  %||% 0) * violence +
          
          # interactions with sex (X = 1)
          (coefs[["sex:H1GI1Y"]]    %||% 0) * 1 * H1GI1Y +
          (coefs[["sex:SES"]]       %||% 0) * 1 * SES +
          (coefs[["sex:happyP"]]    %||% 0) * 1 * happyP +
          (coefs[["sex:drinkP_bin"]]%||% 0) * 1 * drinkP_bin +
          (coefs[["sex:health_bin"]]%||% 0) * 1 * health_bin +
          (coefs[["sex:violence"]]  %||% 0) * 1 * violence +
          
          # covariate–covariate interactions
          (coefs[["H1GI1Y:SES"]]         %||% 0) * H1GI1Y * SES +
          (coefs[["H1GI1Y:happyP"]]      %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]  %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]  %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]    %||% 0) * H1GI1Y * violence +
          (coefs[["SES:happyP"]]         %||% 0) * SES * happyP +
          (coefs[["SES:drinkP_bin"]]     %||% 0) * SES * drinkP_bin +
          (coefs[["SES:health_bin"]]     %||% 0) * SES * health_bin +
          (coefs[["SES:violence"]]       %||% 0) * SES * violence +
          (coefs[["happyP:drinkP_bin"]]  %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]  %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]    %||% 0) * happyP * violence +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
          # residual
          rmse * UY
      )
  })
  
  #----------------------------------------------------------------------
  # Step 3: Model for the Total Causal Effect (TCE), Y ~ X
  #----------------------------------------------------------------------
  
  model <- lm(
    as.formula(
      paste0(
        outcome_var, " ~ ",
        "sex + ",
        "H1GI1Y + SES + ",
        "sex:H1GI1Y + ",
        "sex:SES + ",
        "H1GI1Y:SES"
      )
    ),
    data = dat
  )
  
  coefs <- coef(model)
  
  #----------------------------------------------------------------------
  # Step 3A: Variables for Dep (Y), when X=0 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_0 <- sym(paste0(outcome_var, "_0")) # when X=0
  dat_expanded <- dat_expanded %>%
    mutate(
      !!output_var_0 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[["sex"]] %||% 0) * 0 +
        (coefs[["H1GI1Y"]]    %||% 0) * H1GI1Y +
        (coefs[["SES"]]       %||% 0) * SES +
        (coefs[["sex:H1GI1Y"]]    %||% 0) * H1GI1Y * 0 +
        (coefs[["sex:SES"]]       %||% 0) * SES     * 0 +
        (coefs[["H1GI1Y:SES"]]       %||% 0) * H1GI1Y * SES
    )
  
  #----------------------------------------------------------------------
  # Step 3B: Variables for Dep (Y), when X=1 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_1 <- sym(paste0(outcome_var, "_1")) # when X=1
  dat_expanded <- dat_expanded %>%
    mutate(
      !!output_var_1 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[["sex"]] %||% 0) +
        (coefs[["H1GI1Y"]]    %||% 0) * H1GI1Y +
        (coefs[["SES"]]       %||% 0) * SES +
        (coefs[["sex:H1GI1Y"]]    %||% 0) * 1 * H1GI1Y +
        (coefs[["sex:SES"]]       %||% 0) * 1 * SES +
        (coefs[["H1GI1Y:SES"]]       %||% 0) * H1GI1Y * SES
    )
  
  #--------------------------------------------
  # Step 4: Compute differences: me, tce, dif
  #--------------------------------------------
  
  walk(mediators, function(m) {
    me_var <- sym(paste0("me_", "sex_", m))
    tce_var <- sym(paste0("tce_", "sex_", m))
    dif_var <- sym(paste0("dif_", "sex_", m))
    dep_00 <- sym(paste0(outcome_var, "_00_", m))
    dep_10 <- sym(paste0(outcome_var, "_10_", m))
    dep_0  <- sym(paste0(outcome_var, "_0"))
    dep_1  <- sym(paste0(outcome_var, "_1"))
    
    dat_expanded <<- dat_expanded %>% # calculated for each expanded dataset
      mutate(
        !!me_var := !!dep_10 - !!dep_00,
        !!tce_var := !!dep_1 - !!dep_0,
        !!dif_var := !!tce_var - !!me_var
      )
  })
  
  #-----------------------------------------------------------------
  # Calulate mean me, tce, and dif across all expanded data sets
  #-----------------------------------------------------------------

  effects <- map(mediators, function(m) {
    c(
      ME = mean(dat_expanded[[paste0("me_sex_", m)]], na.rm = TRUE),
      TCE = mean(dat_expanded[[paste0("tce_sex_", m)]], na.rm = TRUE),
      DIF = mean(dat_expanded[[paste0("dif_sex_", m)]], na.rm = TRUE)
    )
  })
  
  names_list <- unlist(map2(mediators, seq_along(mediators),
                            ~ c(paste0("ME", .y), paste0("TCE", .y), paste0("DIF", .y))))
  out <- unlist(effects)
  names(out) <- names_list
  return(out)
}

#----------------------------------
# Run Bootstrapping (outer loop)
#----------------------------------

for (wave in names(waves)) {
  cat("\nBootstrapping wave:", wave, "\n")

  # This re-samples individuals 1000 times to capture sampling variability and calculate confidence intervals 
  boot_results <- boot(
    data = waves[[wave]]$data,
    statistic = function(d, i) bootstrap_stat(d, i, wave_name = wave, mediators = mediators), # runs Monte Carlo expansions within each bootstrap
    R = 1000
  )

  estimates <- boot_results$t
  col_names <- names(boot_results$t0)
  
  # Calculate mean, SE, and percentile confidence intervals
  results <- tibble(
    Effect  = col_names,
    Mean    = boot_results$t0,
    SE      = apply(estimates, 2, sd, na.rm = TRUE),
    CI_Low  = apply(estimates, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_High = apply(estimates, 2, quantile, probs = 0.975, na.rm = TRUE),
    Wave    = wave
  ) %>%
    # Extract mediator number and effect type (ME/TCE/DIF)
    mutate(
      Mediator_Num = as.numeric(gsub("[^0-9]", "", Effect)),
      Effect_Type  = gsub("[0-9]", "", Effect),
      Mediator     = mediators[Mediator_Num],
      Wave         = wave
    ) %>%
    dplyr::select(Wave, Mediator, Effect_Type, Mean, SE, CI_Low, CI_High) # final table

  saveRDS(results, file = paste0("/home/thom1336/healthdis/data/sensitivity/race/confound/sex/bootstrap_results_", wave, race, ".rds"))
  cat("Saved results for wave:", wave, "\n")
  }
}
