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

# list of ses variables
ses_vars <- c("SESlow", "SESmid") 

# create sex specific datasets
dat_pheno_W2_male <- dat_pheno_W2 %>% dplyr::filter(sex == 0)
dat_pheno_W2_female <- dat_pheno_W2 %>% dplyr::filter(sex == 1)

dat_pheno_W4_male <- dat_pheno_W4 %>% dplyr::filter(sex == 0)
dat_pheno_W4_female <- dat_pheno_W4 %>% dplyr::filter(sex == 1)

######----------------------- Runs a loop for each sex-specific dataset -----------------------#######

for(sex in c("male", "female")){
  
  # list of waves and their corresponding data and outcome names
  waves <- list(
    W2 = list(data = get(paste0("dat_pheno_W2_", sex)), outcome = "CESD_W2"),
    W4 = list(data = get(paste0("dat_pheno_W4_", sex)), outcome = "CESD_W4")
  )

# function to run counterfactual calculations and Monte Carlo simulations within bootstrapping 
bootstrap_stat <- function(data, indices, wave_name, mediators, ses_var) {
  
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
  model_list_reg <- list()
  
  #----------------------------------------------------
  # Step 1: Variable for social (M) when ses_var = 0
  #----------------------------------------------------
  
  for (m in mediators) {
    # fit model
    modelM <- lm(as.formula(paste(m, "~", ses_var)), data = dat)  # fit on original bootstrapped sample
    
    # store model in list
    model_list_reg[[paste0("model_", m, "_", ses_var)]] <- modelM
    
    # extract intercept and residual SD
    intercept <- coef(modelM)[["(Intercept)"]]
    resid_sd <- summary(modelM)$sigma
    
    # simulate mediator under ses_var = 0 (on expanded)
    new_var <- paste0(m, "_", ses_var, "0")
    dat_expanded[[new_var]] <- intercept + resid_sd * dat_expanded$UM # applied to all expanded simulations
    
  }
  
  #------------------------------------
  # Step 2: Model for Dep (Y)
  #------------------------------------
  
  walk(mediators, function(m) {
    # construct formula for the model
    formula_str <- paste0(
      outcome_var, " ~ ",                                           # y
      ses_var, " + ",                                               # x
      m, " + ",                                                     # m
      m, ":", ses_var, " + ",                                       # m*x
      m, "2 + ",                                                    # m2
      
      # covariates (main effects)
      "H1GI1Y + ",                                                  # c1
      "black_cov + natam_cov + asian_cov + ",                       # c2 (dummies)
      "happyP + drinkP_bin + health_bin + violence + ",             # new confounders
      
      # ses_var × covariates
      ses_var, ":H1GI1Y + ",
      ses_var, ":black_cov + ",
      ses_var, ":natam_cov + ",
      ses_var, ":asian_cov + ",
      ses_var, ":happyP + ",
      ses_var, ":drinkP_bin + ",
      ses_var, ":health_bin + ",
      ses_var, ":violence + ",
      
      # covariate–covariate interactions
      "H1GI1Y:black_cov + H1GI1Y:natam_cov + H1GI1Y:asian_cov + ",
      "H1GI1Y:happyP + ",
      "H1GI1Y:drinkP_bin + ",
      "H1GI1Y:health_bin + ",
      "H1GI1Y:violence + ",
    #  "happyP:black_cov + happyP:natam_cov + happyP:asian_cov + ",
      "happyP:drinkP_bin + ",
      "happyP:health_bin + ",
      "happyP:violence + ",
    #  "drinkP_bin:black_cov + drinkP_bin:natam_cov + drinkP_bin:asian_cov + ",
      "drinkP_bin:health_bin + ",
      "drinkP_bin:violence + ",
    #  "health_bin:black_cov + health_bin:natam_cov + health_bin:asian_cov + ",
      "health_bin:violence"
    #  "violence:black_cov + violence:natam_cov + violence:asian_cov"
    )
    
    formula_obj <- as.formula(formula_str)
    modelY <- lm(formula_obj, data = dat)  # fit outcome model on original sample
    model_nameY <- paste0("model_", outcome_var, "_", ses_var, "_", m)
    assign(model_nameY, modelY, envir = .GlobalEnv)
    
  })
  
  #-----------------------------------------------------------------------
  # Step 2A: Variable for Y, when ses_var = 0 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", ses_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    ses0_var <- sym(paste0(m, "_", ses_var, "0"))
    output_var_00 <- sym(paste0(outcome_var, "_00_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_00 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # X = ses_var = 0
          (coefs[[ses_var]] %||% 0) * 0 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                          * !!ses0_var +
          (coefs[[paste0(ses_var, ":", m)]] %||% 0)    * 0 * !!ses0_var +
          (coefs[[paste0(m, "2")]] %||% 0)             * (!!ses0_var)^2 +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
          (coefs[["black_cov"]]  %||% 0) * black_cov +
          (coefs[["natam_cov"]]  %||% 0) * natam_cov +
          (coefs[["asian_cov"]]  %||% 0) * asian_cov +
          (coefs[["happyP"]]     %||% 0) * happyP +
          (coefs[["drinkP_bin"]] %||% 0) * drinkP_bin +
          (coefs[["health_bin"]] %||% 0) * health_bin +
          (coefs[["violence"]]   %||% 0) * violence +
          
          # ses_var × covariate interactions (X = 0 → vanish)
          (coefs[[paste0(ses_var, ":H1GI1Y")]]     %||% 0) * 0 * H1GI1Y +
          (coefs[[paste0(ses_var, ":black_cov")]]  %||% 0) * 0 * black_cov +
          (coefs[[paste0(ses_var, ":natam_cov")]]  %||% 0) * 0 * natam_cov +
          (coefs[[paste0(ses_var, ":asian_cov")]]  %||% 0) * 0 * asian_cov +
          (coefs[[paste0(ses_var, ":happyP")]]     %||% 0) * 0 * happyP +
          (coefs[[paste0(ses_var, ":drinkP_bin")]] %||% 0) * 0 * drinkP_bin +
          (coefs[[paste0(ses_var, ":health_bin")]] %||% 0) * 0 * health_bin +
          (coefs[[paste0(ses_var, ":violence")]]   %||% 0) * 0 * violence +
          
          # covariate–covariate interactions (expanded)
          (coefs[["H1GI1Y:black_cov"]]    %||% 0) * H1GI1Y * black_cov +
          (coefs[["H1GI1Y:natam_cov"]]    %||% 0) * H1GI1Y * natam_cov +
          (coefs[["H1GI1Y:asian_cov"]]    %||% 0) * H1GI1Y * asian_cov +
          (coefs[["H1GI1Y:happyP"]]       %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]   %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]   %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]     %||% 0) * H1GI1Y * violence +
          
        #  (coefs[["happyP:black_cov"]]    %||% 0) * happyP * black_cov +
        #  (coefs[["happyP:natam_cov"]]    %||% 0) * happyP * natam_cov +
        #  (coefs[["happyP:asian_cov"]]    %||% 0) * happyP * asian_cov +
          (coefs[["happyP:drinkP_bin"]]   %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]   %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]     %||% 0) * happyP * violence +
          
         #  (coefs[["drinkP_bin:black_cov"]]  %||% 0) * drinkP_bin * black_cov +
         #  (coefs[["drinkP_bin:natam_cov"]]  %||% 0) * drinkP_bin * natam_cov +
         #  (coefs[["drinkP_bin:asian_cov"]]  %||% 0) * drinkP_bin * asian_cov +
           (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
           (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          
       #   (coefs[["health_bin:black_cov"]]  %||% 0) * health_bin * black_cov +
       #   (coefs[["health_bin:natam_cov"]]  %||% 0) * health_bin * natam_cov +
       #   (coefs[["health_bin:asian_cov"]]  %||% 0) * health_bin * asian_cov +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
       #   (coefs[["violence:black_cov"]]    %||% 0) * violence * black_cov +
       #   (coefs[["violence:natam_cov"]]    %||% 0) * violence * natam_cov +
       #   (coefs[["violence:asian_cov"]]    %||% 0) * violence * asian_cov +
          
          # residual
          rmse * UY
      )
  })
  
  #----------------------------------------------------------------------
  # Step 2B: Variable for Y, when ses_var = 1 and M=m_0, using model from Step 2 
  #----------------------------------------------------------------------
  
  walk(mediators, function(m) {
    model_obj <- get(paste0("model_", outcome_var, "_", ses_var, "_", m))
    coefs <- coef(model_obj)
    rmse <- summary(model_obj)$sigma
    
    ses0_var   <- sym(paste0(m, "_", ses_var, "0"))
    output_var_10 <- sym(paste0(outcome_var, "_10_", m))
    
    dat_expanded <<- dat_expanded %>%
      mutate(
        !!output_var_10 :=
          # intercept
          (coefs[["(Intercept)"]] %||% 0) +
          
          # X = ses_var = 1
          (coefs[[ses_var]] %||% 0) * 1 +
          
          # mediator terms
          (coefs[[m]] %||% 0)                          * !!ses0_var +
          (coefs[[paste0(ses_var, ":", m)]] %||% 0)    * 1 * !!ses0_var +
          (coefs[[paste0(m, "2")]] %||% 0)             * (!!ses0_var)^2 +
          
          # covariate main effects
          (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
          (coefs[["black_cov"]]  %||% 0) * black_cov +
          (coefs[["natam_cov"]]  %||% 0) * natam_cov +
          (coefs[["asian_cov"]]  %||% 0) * asian_cov +
          (coefs[["happyP"]]     %||% 0) * happyP +
          (coefs[["drinkP_bin"]] %||% 0) * drinkP_bin +
          (coefs[["health_bin"]] %||% 0) * health_bin +
          (coefs[["violence"]]   %||% 0) * violence +
          
          # ses_var × covariate interactions (X = 1)
          (coefs[[paste0(ses_var, ":H1GI1Y")]]     %||% 0) * 1 * H1GI1Y +
          (coefs[[paste0(ses_var, ":black_cov")]]  %||% 0) * 1 * black_cov +
          (coefs[[paste0(ses_var, ":natam_cov")]]  %||% 0) * 1 * natam_cov +
          (coefs[[paste0(ses_var, ":asian_cov")]]  %||% 0) * 1 * asian_cov +
          (coefs[[paste0(ses_var, ":happyP")]]     %||% 0) * 1 * happyP +
          (coefs[[paste0(ses_var, ":drinkP_bin")]] %||% 0) * 1 * drinkP_bin +
          (coefs[[paste0(ses_var, ":health_bin")]] %||% 0) * 1 * health_bin +
          (coefs[[paste0(ses_var, ":violence")]]   %||% 0) * 1 * violence +
          
          # covariate–covariate interactions (expanded)
          (coefs[["H1GI1Y:black_cov"]]    %||% 0) * H1GI1Y * black_cov +
          (coefs[["H1GI1Y:natam_cov"]]    %||% 0) * H1GI1Y * natam_cov +
          (coefs[["H1GI1Y:asian_cov"]]    %||% 0) * H1GI1Y * asian_cov +
          (coefs[["H1GI1Y:happyP"]]       %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]   %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]   %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]     %||% 0) * H1GI1Y * violence +
          
       #   (coefs[["happyP:black_cov"]]    %||% 0) * happyP * black_cov +
       #   (coefs[["happyP:natam_cov"]]    %||% 0) * happyP * natam_cov +
       #   (coefs[["happyP:asian_cov"]]    %||% 0) * happyP * asian_cov +
          (coefs[["happyP:drinkP_bin"]]   %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]   %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]     %||% 0) * happyP * violence +
          
      #    (coefs[["drinkP_bin:black_cov"]]  %||% 0) * drinkP_bin * black_cov +
      #    (coefs[["drinkP_bin:natam_cov"]]  %||% 0) * drinkP_bin * natam_cov +
      #    (coefs[["drinkP_bin:asian_cov"]]  %||% 0) * drinkP_bin * asian_cov +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          
       #   (coefs[["health_bin:black_cov"]]  %||% 0) * health_bin * black_cov +
       #   (coefs[["health_bin:natam_cov"]]  %||% 0) * health_bin * natam_cov +
       #   (coefs[["health_bin:asian_cov"]]  %||% 0) * health_bin * asian_cov +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          
       #   (coefs[["violence:black_cov"]]    %||% 0) * violence * black_cov +
       #   (coefs[["violence:natam_cov"]]    %||% 0) * violence * natam_cov +
       #   (coefs[["violence:asian_cov"]]    %||% 0) * violence * asian_cov +
          
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
        ses_var, " + ",
        "H1GI1Y + black_cov + natam_cov + asian_cov + ",
        ses_var, ":H1GI1Y + ",
        ses_var, ":black_cov + ",
        ses_var, ":natam_cov + ",
        ses_var, ":asian_cov + ",
        "H1GI1Y:black_cov + H1GI1Y:natam_cov + H1GI1Y:asian_cov"
      )
    ),
    data = dat
  ) # computed on original data
  
  coefs <- coef(model)
  
  #----------------------------------------------------------------------
  # Step 3A: Variables for Dep (Y), when ses_var = 0 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_0 <- sym(paste0(outcome_var, "_0")) # when ses_var = 0
  dat_expanded <- dat_expanded %>%
    mutate(
      !!output_var_0 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[[ses_var]] %||% 0) * 0 +
        (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
        (coefs[["black_cov"]]  %||% 0) * black_cov +
        (coefs[["natam_cov"]]  %||% 0) * natam_cov +
        (coefs[["asian_cov"]]  %||% 0) * asian_cov +
        (coefs[[paste0(ses_var, ":H1GI1Y")]]  %||% 0) * 0 * H1GI1Y +
        (coefs[[paste0(ses_var, ":black_cov")]] %||% 0) * 0 * black_cov +
        (coefs[[paste0(ses_var, ":natam_cov")]] %||% 0) * 0 * natam_cov +
        (coefs[[paste0(ses_var, ":asian_cov")]] %||% 0) * 0 * asian_cov +
        (coefs[["H1GI1Y:black_cov"]] %||% 0) * H1GI1Y * black_cov +
        (coefs[["H1GI1Y:natam_cov"]] %||% 0) * H1GI1Y * natam_cov +
        (coefs[["H1GI1Y:asian_cov"]] %||% 0) * H1GI1Y * asian_cov 
    )
  
  #----------------------------------------------------------------------
  # Step 3B: Variables for Dep (Y), when ses_var = 1 (no mediator involved) 
  #----------------------------------------------------------------------
  
  output_var_1 <- sym(paste0(outcome_var, "_1")) # when ses_var = 1
  dat_expanded <- dat_expanded %>%
    mutate(
      !!output_var_1 :=
        (coefs[["(Intercept)"]] %||% 0) +
        (coefs[[ses_var]] %||% 0) +
        (coefs[["H1GI1Y"]]     %||% 0) * H1GI1Y +
        (coefs[["black_cov"]]  %||% 0) * black_cov +
        (coefs[["natam_cov"]]  %||% 0) * natam_cov +
        (coefs[["asian_cov"]]  %||% 0) * asian_cov +
        (coefs[[paste0(ses_var, ":H1GI1Y")]]  %||% 0) * 1 * H1GI1Y +
        (coefs[[paste0(ses_var, ":black_cov")]] %||% 0) * 1 * black_cov +
        (coefs[[paste0(ses_var, ":natam_cov")]] %||% 0) * 1 * natam_cov +
        (coefs[[paste0(ses_var, ":asian_cov")]] %||% 0) * 1 * asian_cov +
        (coefs[["H1GI1Y:black_cov"]] %||% 0) * H1GI1Y * black_cov +
        (coefs[["H1GI1Y:natam_cov"]] %||% 0) * H1GI1Y * natam_cov +
        (coefs[["H1GI1Y:asian_cov"]] %||% 0) * H1GI1Y * asian_cov 
    )
  
  #--------------------------------------------
  # Step 4: Compute differences: de, tce, dif
  #--------------------------------------------
  
  walk(mediators, function(m) {
    me_var <- sym(paste0("me_", ses_var, "_", m))
    tce_var <- sym(paste0("tce_", ses_var, "_", m))
    dif_var <- sym(paste0("dif_", ses_var, "_", m))
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
      ME = mean(dat_expanded[[paste0("me_", ses_var, "_", m)]], na.rm = TRUE),
      TCE = mean(dat_expanded[[paste0("tce_", ses_var, "_", m)]], na.rm = TRUE),
      DIF = mean(dat_expanded[[paste0("dif_", ses_var, "_", m)]], na.rm = TRUE)
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
  for (ses_var in ses_vars) {  # loop over each ses variable
    cat("\nBootstrapping wave:", wave, "for ses variable:", ses_var, "\n")
    
    # This re-samples individuals 1000 times to capture sampling variability and calculate confidence intervals 
    boot_results <- boot(
      data = waves[[wave]]$data,
      statistic = function(d, i) bootstrap_stat(d, i, wave_name = wave, mediators = mediators, ses_var = ses_var), # loops through ses_var
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
      Wave    = wave,
      ses    = ses_var  # Add the ses variable to the results
    ) %>%
      # Extract mediator number and effect type (ME/TCE/DIF)
      mutate(
        Mediator_Num = as.numeric(gsub("[^0-9]", "", Effect)),
        Effect_Type  = gsub("[0-9]", "", Effect),
        Mediator     = mediators[Mediator_Num],
        Wave         = wave
      ) %>%
      dplyr::select(Wave, ses, Mediator, Effect_Type, Mean, SE, CI_Low, CI_High) # final table
    
    saveRDS(results, file = paste0("/home/thom1336/healthdis/data/sensitivity/sex/confound/ses/bootstrap_results_", wave, sex, "_", ses_var, ".rds"))
    cat("Saved results for wave:", wave, "and ses variable:", ses_var, "\n")
    }
  }
}
