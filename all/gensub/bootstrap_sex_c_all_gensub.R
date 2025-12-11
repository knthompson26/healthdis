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
load("/home/thom1336/healthdis/data/dat_pheno_pgs_W2.RData")
load("/home/thom1336/healthdis/data/dat_pheno_pgs_W4.RData")

# list of mediators
mediators <- c("family", "school", "friends", "neigh", "religion", "social")

# list of waves and their corresponding data and outcome names
waves <- list(
  W2 = list(data = dat_pheno_pgs_W2, outcome = "CESD_W2"),
  W4 = list(data = dat_pheno_pgs_W4, outcome = "CESD_W4")
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
      "black_cov + natam_cov + asian_cov + ",                      # c3 (dummies)
      "happyP + drinkP_bin + health_bin + violence + ",            # parent nd health
      "pgs + ",                                                    # new confounder - pgs
      
      # sex × covariates
      "sex:H1GI1Y + ",                                             # c1*x
      "sex:SES + ",                                                # c2*x 
      "sex:black_cov + sex:natam_cov + sex:asian_cov + ",          # c3*x 
      "sex:happyP + ",
      "sex:drinkP_bin + ",
      "sex:health_bin + ",
      "sex:violence + ",
      "sex:pgs + ",
      
      # coavriate x covariate
      "H1GI1Y:SES + ",                                             
      "H1GI1Y:black_cov + H1GI1Y:natam_cov + H1GI1Y:asian_cov + ", 
      "H1GI1Y:happyP + ",
      "H1GI1Y:drinkP_bin + ",
      "H1GI1Y:health_bin + ",
      "H1GI1Y:violence + ",
      "H1GI1Y:pgs + ",
      "SES:black_cov + SES:natam_cov + SES:asian_cov +",
      "SES:happyP + ",
      "SES:drinkP_bin + ",
      "SES:health_bin + ",
      "SES:violence + ",
      "SES:pgs + ",
   #   "happyP:black_cov + happyP:natam_cov + happyP:asian_cov + ",
      "happyP:drinkP_bin + ",
      "happyP:health_bin + ",
      "happyP:violence + ",
      "happyP:pgs + ",
   #   "drinkP_bin:black_cov + drinkP_bin:natam_cov + drinkP_bin:asian_cov + ",
      "drinkP_bin:health_bin + ",
      "drinkP_bin:violence + ",
      "drinkP_bin:pgs + ",
   #   "health_bin:black_cov + health_bin:natam_cov + health_bin:asian_cov + ",
      "health_bin:violence +",
      "health_bin:pgs + ",
  #    "violence:black_cov + violence:natam_cov + violence:asian_cov + ",
      "violence:pgs"
   #   "pgs:black_cov + pgs:natam_cov + pgs:asian_cov"
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
          (coefs[["black_cov"]] %||% 0) * black_cov +
          (coefs[["natam_cov"]] %||% 0) * natam_cov +
          (coefs[["asian_cov"]] %||% 0) * asian_cov +
          (coefs[["happyP"]]    %||% 0) * happyP +
          (coefs[["drinkP_bin"]]%||% 0) * drinkP_bin +
          (coefs[["health_bin"]]%||% 0) * health_bin +
          (coefs[["violence"]]  %||% 0) * violence + 
          (coefs[["pgs"]]  %||% 0) * pgs + 
          
          # interactions with sex (X = 0 → all vanish)
          (coefs[["sex:H1GI1Y"]]    %||% 0) * H1GI1Y   * 0 +
          (coefs[["sex:SES"]]       %||% 0) * SES       * 0 +
          (coefs[["sex:black_cov"]] %||% 0) * black_cov * 0 +
          (coefs[["sex:natam_cov"]] %||% 0) * natam_cov * 0 +
          (coefs[["sex:asian_cov"]] %||% 0) * asian_cov * 0 +
          (coefs[["sex:happyP"]]    %||% 0) * happyP    * 0 +
          (coefs[["sex:drinkP_bin"]]%||% 0) * drinkP_bin * 0 +
          (coefs[["sex:health_bin"]]%||% 0) * health_bin * 0 +
          (coefs[["sex:violence"]]  %||% 0) * violence  * 0 + 
          (coefs[["sex:pgs"]]  %||% 0) * pgs  * 0 + 
          
          # covariate–covariate interactions
          (coefs[["H1GI1Y:SES"]]         %||% 0) * H1GI1Y * SES +
          (coefs[["H1GI1Y:black_cov"]]   %||% 0) * H1GI1Y * black_cov +
          (coefs[["H1GI1Y:natam_cov"]]   %||% 0) * H1GI1Y * natam_cov +
          (coefs[["H1GI1Y:asian_cov"]]   %||% 0) * H1GI1Y * asian_cov +
          (coefs[["H1GI1Y:happyP"]]      %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]  %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]  %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]    %||% 0) * H1GI1Y * violence +
          (coefs[["H1GI1Y:pgs"]]    %||% 0) * H1GI1Y * pgs +
          
          (coefs[["SES:black_cov"]]      %||% 0) * SES * black_cov +
          (coefs[["SES:natam_cov"]]      %||% 0) * SES * natam_cov +
          (coefs[["SES:asian_cov"]]      %||% 0) * SES * asian_cov +
          (coefs[["SES:happyP"]]         %||% 0) * SES * happyP +
          (coefs[["SES:drinkP_bin"]]     %||% 0) * SES * drinkP_bin +
          (coefs[["SES:health_bin"]]     %||% 0) * SES * health_bin +
          (coefs[["SES:violence"]]       %||% 0) * SES * violence +
          (coefs[["SES:pgs"]]       %||% 0) * SES * pgs +
          
     #     (coefs[["happyP:black_cov"]]   %||% 0) * happyP * black_cov +
     #     (coefs[["happyP:natam_cov"]]   %||% 0) * happyP * natam_cov +
     #     (coefs[["happyP:asian_cov"]]   %||% 0) * happyP * asian_cov +
          (coefs[["happyP:drinkP_bin"]]  %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]  %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]    %||% 0) * happyP * violence +
          (coefs[["happyP:pgs"]]    %||% 0) * happyP * pgs +
          
     #     (coefs[["drinkP_bin:black_cov"]]  %||% 0) * drinkP_bin * black_cov +
     #     (coefs[["drinkP_bin:natam_cov"]]  %||% 0) * drinkP_bin * natam_cov +
     #     (coefs[["drinkP_bin:asian_cov"]]  %||% 0) * drinkP_bin * asian_cov +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          (coefs[["drinkP_bin:pgs"]]   %||% 0) * drinkP_bin * pgs +
          
      #    (coefs[["health_bin:black_cov"]]  %||% 0) * health_bin * black_cov +
      #    (coefs[["health_bin:natam_cov"]]  %||% 0) * health_bin * natam_cov +
      #    (coefs[["health_bin:asian_cov"]]  %||% 0) * health_bin * asian_cov +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          (coefs[["health_bin:pgs"]]   %||% 0) * health_bin * pgs +
          
      #    (coefs[["violence:black_cov"]]    %||% 0) * violence * black_cov +
      #    (coefs[["violence:natam_cov"]]    %||% 0) * violence * natam_cov +
      #    (coefs[["violence:asian_cov"]]    %||% 0) * violence * asian_cov +
          (coefs[["violence:pgs"]]    %||% 0) * violence * pgs +
          
      #    (coefs[["pgs:black_cov"]]    %||% 0) * pgs * black_cov +
      #    (coefs[["pgs:natam_cov"]]    %||% 0) * pgs * natam_cov +
      #    (coefs[["pgs:asian_cov"]]    %||% 0) * pgs * asian_cov +
          
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
          (coefs[["black_cov"]] %||% 0) * black_cov +
          (coefs[["natam_cov"]] %||% 0) * natam_cov +
          (coefs[["asian_cov"]] %||% 0) * asian_cov +
          (coefs[["happyP"]]    %||% 0) * happyP +
          (coefs[["drinkP_bin"]]%||% 0) * drinkP_bin +
          (coefs[["health_bin"]]%||% 0) * health_bin +
          (coefs[["violence"]]  %||% 0) * violence +
          (coefs[["pgs"]]  %||% 0) * pgs +
          
          # interactions with sex (X = 1)
          (coefs[["sex:H1GI1Y"]]    %||% 0) * 1 * H1GI1Y +
          (coefs[["sex:SES"]]       %||% 0) * 1 * SES +
          (coefs[["sex:black_cov"]] %||% 0) * 1 * black_cov +
          (coefs[["sex:natam_cov"]] %||% 0) * 1 * natam_cov +
          (coefs[["sex:asian_cov"]] %||% 0) * 1 * asian_cov +
          (coefs[["sex:happyP"]]    %||% 0) * 1 * happyP +
          (coefs[["sex:drinkP_bin"]]%||% 0) * 1 * drinkP_bin +
          (coefs[["sex:health_bin"]]%||% 0) * 1 * health_bin +
          (coefs[["sex:violence"]]  %||% 0) * 1 * violence +
          (coefs[["sex:pgs"]]  %||% 0) * 1 * pgs +
          
          # covariate–covariate interactions
          (coefs[["H1GI1Y:SES"]]         %||% 0) * H1GI1Y * SES +
          (coefs[["H1GI1Y:black_cov"]]   %||% 0) * H1GI1Y * black_cov +
          (coefs[["H1GI1Y:natam_cov"]]   %||% 0) * H1GI1Y * natam_cov +
          (coefs[["H1GI1Y:asian_cov"]]   %||% 0) * H1GI1Y * asian_cov +
          (coefs[["H1GI1Y:happyP"]]      %||% 0) * H1GI1Y * happyP +
          (coefs[["H1GI1Y:drinkP_bin"]]  %||% 0) * H1GI1Y * drinkP_bin +
          (coefs[["H1GI1Y:health_bin"]]  %||% 0) * H1GI1Y * health_bin +
          (coefs[["H1GI1Y:violence"]]    %||% 0) * H1GI1Y * violence +
          (coefs[["H1GI1Y:pgs"]]    %||% 0) * H1GI1Y * pgs +
          
          (coefs[["SES:black_cov"]]      %||% 0) * SES * black_cov +
          (coefs[["SES:natam_cov"]]      %||% 0) * SES * natam_cov +
          (coefs[["SES:asian_cov"]]      %||% 0) * SES * asian_cov +
          (coefs[["SES:happyP"]]         %||% 0) * SES * happyP +
          (coefs[["SES:drinkP_bin"]]     %||% 0) * SES * drinkP_bin +
          (coefs[["SES:health_bin"]]     %||% 0) * SES * health_bin +
          (coefs[["SES:violence"]]       %||% 0) * SES * violence +
          (coefs[["SES:pgs"]]       %||% 0) * SES * pgs +
          
     #     (coefs[["happyP:black_cov"]]   %||% 0) * happyP * black_cov +
     #     (coefs[["happyP:natam_cov"]]   %||% 0) * happyP * natam_cov +
     #     (coefs[["happyP:asian_cov"]]   %||% 0) * happyP * asian_cov +
          (coefs[["happyP:drinkP_bin"]]  %||% 0) * happyP * drinkP_bin +
          (coefs[["happyP:health_bin"]]  %||% 0) * happyP * health_bin +
          (coefs[["happyP:violence"]]    %||% 0) * happyP * violence +
          (coefs[["happyP:pgs"]]    %||% 0) * happyP * pgs +
          
    #      (coefs[["drinkP_bin:black_cov"]]  %||% 0) * drinkP_bin * black_cov +
    #      (coefs[["drinkP_bin:natam_cov"]]  %||% 0) * drinkP_bin * natam_cov +
    #      (coefs[["drinkP_bin:asian_cov"]]  %||% 0) * drinkP_bin * asian_cov +
          (coefs[["drinkP_bin:health_bin"]] %||% 0) * drinkP_bin * health_bin +
          (coefs[["drinkP_bin:violence"]]   %||% 0) * drinkP_bin * violence +
          (coefs[["drinkP_bin:pgs"]]   %||% 0) * drinkP_bin * pgs +
          
     #     (coefs[["health_bin:black_cov"]]  %||% 0) * health_bin * black_cov +
     #     (coefs[["health_bin:natam_cov"]]  %||% 0) * health_bin * natam_cov +
     #     (coefs[["health_bin:asian_cov"]]  %||% 0) * health_bin * asian_cov +
          (coefs[["health_bin:violence"]]   %||% 0) * health_bin * violence +
          (coefs[["health_bin:pgs"]]   %||% 0) * health_bin * pgs +
          
     #     (coefs[["violence:black_cov"]]    %||% 0) * violence * black_cov +
     #     (coefs[["violence:natam_cov"]]    %||% 0) * violence * natam_cov +
     #     (coefs[["violence:asian_cov"]]    %||% 0) * violence * asian_cov +
          (coefs[["violence:pgs"]]    %||% 0) * violence * pgs +
          
     #     (coefs[["pgs:black_cov"]]    %||% 0) * pgs * black_cov +
     #     (coefs[["pgs:natam_cov"]]    %||% 0) * pgs * natam_cov +
     #     (coefs[["pgs:asian_cov"]]    %||% 0) * pgs * asian_cov +
          
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
        "H1GI1Y + SES + black_cov + natam_cov + asian_cov + pgs +",
        "sex:H1GI1Y + sex:SES + sex:black_cov + sex:natam_cov + sex:asian_cov + sex:pgs +",
        "H1GI1Y:SES + H1GI1Y:black_cov + H1GI1Y:natam_cov + H1GI1Y:asian_cov + H1GI1Y:pgs +",
        "SES:black_cov + SES:natam_cov + SES:asian_cov + SES:pgs"
       # "pgs:black_cov + pgs:natam_cov + pgs:asian_cov"
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
        (coefs[["black_cov"]] %||% 0) * black_cov +
        (coefs[["natam_cov"]] %||% 0) * natam_cov +
        (coefs[["asian_cov"]] %||% 0) * asian_cov +
        (coefs[["pgs"]] %||% 0) * pgs +
        (coefs[["sex:H1GI1Y"]]    %||% 0) * H1GI1Y * 0 +
        (coefs[["sex:SES"]]       %||% 0) * SES     * 0 +
        (coefs[["sex:black_cov"]] %||% 0) * black_cov * 0 +
        (coefs[["sex:natam_cov"]] %||% 0) * natam_cov * 0 +
        (coefs[["sex:asian_cov"]] %||% 0) * asian_cov * 0 +
        (coefs[["sex:pgs"]] %||% 0) * pgs * 0 +
        (coefs[["H1GI1Y:SES"]]       %||% 0) * H1GI1Y * SES +
        (coefs[["H1GI1Y:black_cov"]] %||% 0) * H1GI1Y * black_cov +
        (coefs[["H1GI1Y:natam_cov"]] %||% 0) * H1GI1Y * natam_cov +
        (coefs[["H1GI1Y:asian_cov"]] %||% 0) * H1GI1Y * asian_cov +
        (coefs[["H1GI1Y:pgs"]] %||% 0) * H1GI1Y * pgs +
        (coefs[["SES:black_cov"]]    %||% 0) * SES * black_cov +
        (coefs[["SES:natam_cov"]]    %||% 0) * SES * natam_cov +
        (coefs[["SES:asian_cov"]]    %||% 0) * SES * asian_cov +
        (coefs[["SES:pgs"]]    %||% 0) * SES * pgs 
     #   (coefs[["pgs:black_cov"]]    %||% 0) * pgs * black_cov +
     #   (coefs[["pgs:natam_cov"]]    %||% 0) * pgs * natam_cov +
     #   (coefs[["pgs:asian_cov"]]    %||% 0) * pgs * asian_cov
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
        (coefs[["black_cov"]] %||% 0) * black_cov +
        (coefs[["natam_cov"]] %||% 0) * natam_cov +
        (coefs[["asian_cov"]] %||% 0) * asian_cov +
        (coefs[["pgs"]] %||% 0) * pgs +
        (coefs[["sex:H1GI1Y"]]    %||% 0) * 1 * H1GI1Y +
        (coefs[["sex:SES"]]       %||% 0) * 1 * SES +
        (coefs[["sex:black_cov"]] %||% 0) * 1 * black_cov +
        (coefs[["sex:natam_cov"]] %||% 0) * 1 * natam_cov +
        (coefs[["sex:asian_cov"]] %||% 0) * 1 * asian_cov +
        (coefs[["sex:pgs"]] %||% 0) * 1 * pgs +
        (coefs[["H1GI1Y:SES"]]       %||% 0) * H1GI1Y * SES +
        (coefs[["H1GI1Y:black_cov"]] %||% 0) * H1GI1Y * black_cov +
        (coefs[["H1GI1Y:natam_cov"]] %||% 0) * H1GI1Y * natam_cov +
        (coefs[["H1GI1Y:asian_cov"]] %||% 0) * H1GI1Y * asian_cov +
        (coefs[["H1GI1Y:pgs"]] %||% 0) * H1GI1Y * pgs +
        (coefs[["SES:black_cov"]]    %||% 0) * SES * black_cov +
        (coefs[["SES:natam_cov"]]    %||% 0) * SES * natam_cov +
        (coefs[["SES:asian_cov"]]    %||% 0) * SES * asian_cov +
        (coefs[["SES:pgs"]]    %||% 0) * SES * pgs 
     #   (coefs[["pgs:black_cov"]]    %||% 0) * pgs * black_cov +
     #   (coefs[["pgs:natam_cov"]]    %||% 0) * pgs * natam_cov +
     #   (coefs[["pgs:asian_cov"]]    %||% 0) * pgs * asian_cov
    )
  
  #--------------------------------------------
  # Step 4: Compute differences: de, tce, dif
  #--------------------------------------------
  
  walk(mediators, function(m) {
    de_var <- sym(paste0("de_", "sex_", m))
    tce_var <- sym(paste0("tce_", "sex_", m))
    dif_var <- sym(paste0("dif_", "sex_", m))
    dep_00 <- sym(paste0(outcome_var, "_00_", m))
    dep_10 <- sym(paste0(outcome_var, "_10_", m))
    dep_0  <- sym(paste0(outcome_var, "_0"))
    dep_1  <- sym(paste0(outcome_var, "_1"))
    
    dat_expanded <<- dat_expanded %>% # calculated for each expanded dataset
      mutate(
        !!de_var := !!dep_10 - !!dep_00,
        !!tce_var := !!dep_1 - !!dep_0,
        !!dif_var := !!tce_var - !!de_var
      )
  })
  
  #-----------------------------------------------------------------
  # Calulate mean de, tce, and dif across all expanded data sets
  #-----------------------------------------------------------------

  effects <- map(mediators, function(m) {
    c(
      CDM = mean(dat_expanded[[paste0("de_sex_", m)]], na.rm = TRUE),
      TCE = mean(dat_expanded[[paste0("tce_sex_", m)]], na.rm = TRUE),
      DIF = mean(dat_expanded[[paste0("dif_sex_", m)]], na.rm = TRUE)
    )
  })
  
  names_list <- unlist(map2(mediators, seq_along(mediators),
                            ~ c(paste0("CDM", .y), paste0("TCE", .y), paste0("DIF", .y))))
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
    # Extract mediator number and effect type (CDM/TCE/DIF)
    mutate(
      Mediator_Num = as.numeric(gsub("[^0-9]", "", Effect)),
      Effect_Type  = gsub("[0-9]", "", Effect),
      Mediator     = mediators[Mediator_Num],
      Wave         = wave
    ) %>%
    dplyr::select(Wave, Mediator, Effect_Type, Mean, SE, CI_Low, CI_High) # final table

  saveRDS(results, file = paste0("/home/thom1336/healthdis/data/confound/all/gensub/sex/bootstrap_results_", wave, ".rds"))
 # saveRDS(boot_results, file = paste0("/home/thom1336/healthdis/data/confound/all/gensub/sex/boot_object_", wave, ".rds"))

  cat("Saved results for wave:", wave, "\n")
}
