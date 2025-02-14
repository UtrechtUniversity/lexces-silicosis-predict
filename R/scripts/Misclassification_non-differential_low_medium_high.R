#### About this script ####
## Title: Estimation of the impact of outcome misclassification in a
##        pneumoconiosis diagnostic prediction rule
## Author: Javier Mancilla Galindo, Lützen Portengen 
## ORCiD: https://orcid.org/0000-0002-0718-467X
##
## Purpose: This script follows the procedure described by Zawitowski, et al. 
## https://onlinelibrary.wiley.com/doi/10.1002/sim.7260 to estimate the impact 
## of outcome misclassification in a pneumoconiosis diagnostic prediction rule
## https://oem.bmj.com/lookup/doi/10.1136/oem.2006.027904. Available code 
## in Zawitowski's paper was adapted to simulate data compatible with the 
## summary sample data and prediction model equation in the OEM paper.
## A total of 5000 iterations are made to obtain stable mean AUC estimates. 
##
## Note: Functions by Zawitowski, et al. are found elsewhere. This script is
## sourced into the main quarto markdown file (.qmd) after function loading. 
## Refer to README file of repository for instructions on how to use. 

#### START ####

#### Source functions #### 
### This is disabled by default. Use only when running script independently
### from the main qmd file, by deleting the following "#" character:

# source("R/scripts/Zawistowski_misclassification_functions.R")

#### Set parameters ####

# Set seed for reproducibility
set.seed(2024) # Current year

# Number of repetitions
num_repetitions <- 5000

# Sample characteristics:

# source("R/scripts/sample_characteristics_simulation.R")
  
# Misclassification parameters were obtained from a systematic review
# and meta-analysis of the diagnostic accuracy of CXR against HRCT for 
# https://thorax.bmj.com/lookup/doi/10.1136/thorax-2024-BTSabstracts.23
# and associated data in the repository: 
# https://github.com/pjhowlett/da_silic_cxr

# Values are defined in environment after sourcing into qmd file the script: 
# Howlett_CXR_HRCT_dx_performance.R 

# Get beta parameters for the beta distribution of sensitivity
beta_params <- get.beta.par(
  p = c(0.025, 0.975),                      # For 95% CI
  q = c(pooled_ci_Sn[1], pooled_ci_Sn[2]),  # Pooled sensitivity 95%CI
  m = pooled_Sn,                            # Pooled sensitivity
  plot = FALSE,
  show = FALSE
)

alpha <- beta_params[1]
beta <- beta_params[2]

# Initialize vectors to store results
ROC_observed <- vector(length = num_repetitions)
ROC_observed_corrected <- vector(length = num_repetitions)
sum_low_risk <- vector(length = num_repetitions)
sum_medium_risk <- vector(length = num_repetitions)
sum_high_risk <- vector(length = num_repetitions)
sum_outcome_observed <- vector(length = num_repetitions)
sum_outcome_true <- vector(length = num_repetitions)
sum_high_risk_observed <- vector(length = num_repetitions)
sum_high_risk_true <- vector(length = num_repetitions)
sum_medium_risk_observed <- vector(length = num_repetitions)
sum_medium_risk_true <- vector(length = num_repetitions)
sum_low_risk_observed <- vector(length = num_repetitions)
sum_low_risk_true <- vector(length = num_repetitions)
mean_age_low_risk_true <- vector(length = num_repetitions)
mean_age_medium_risk_true <- vector(length = num_repetitions)
mean_age_high_risk_true <- vector(length = num_repetitions)

#### Simulation of data, scores, and AUC estimates for every new sample ####

for (i in 1:num_repetitions) {
  # Show progress every 500 iterations
  if (i %% 500 == 0) {
    progress <- (i/num_repetitions) * 100
    cat(sprintf("Completed %d simulations (%.1f%%)\n", i, progress))
  }
  
  # Sample Sensitivity from Beta distribution from meta-analysis
  Sn_sim <- rbeta(1, alpha, beta)
    
  # Compute Specificity (Sp) using the DOR formula, with corrected DOR
  Sp_sim <- calculate_Sp(Sn_sim, pooled_DOR)
  
  # Misclassification parameters (varying for each iteration)
  g0 <- 1 - Sp_sim  # gamma0 (false-positive rate)
  g1 <- 1 - Sn_sim  # gamma1 (false-negative rate)
  
  # Simulate Data
  age <- rnorm(N, age_mean, age_sd)
  
  X1 <- rbinom(N, 1, prop_1) # Age > 40
  age[X1 == 1] <- pmax(age[X1 == 1], 40) # Ensure age > 40 where X1 = 1
  X2 <- rbinom(N, 1, prop_2) # Current smoker
  X3 <- rbinom(N, 1, prop_3) # High exposed job title
  X4 <- rbinom(N, 1, prop_4) # Work in the construction industry > 15 years
  X5 <- rbinom(N, 1, prop_5) # Feeling unhealthy
  X6 <- rbinom(N, 1, prop_6) # Standardized residual FEV1 <-1.0
  
  # Outcome probability, based on prediction rule equation: 
  p <- exp(
    B0 + B1*X1 + B2*X2 + B3*X3 + B4*X4 + B5*X5 + B6*X6) / 
    (1 + exp(B0 + B1*X1 + B2*X2 + B3*X3 + B4*X4 + B5*X5 + B6*X6)
    ) 
  
  # Simulate observed outcome based on prediction rule formula
  Outcome_observed <- rbinom(N, 1, p)
  
  # Reverse-misclassify outcome (Index test to reference test)
  Outcome_true <- misclassify_reverse(Outcome_observed, Sn_sim, Sp_sim)
  
  # Calculate scores according to the diagnostic prediction rule
  score <- X1*1 + X2*1 + X3*1.5 + X4*1.5 + X5*1.25 + X6*1.25
  
  # Categorize into high and low risk. 
  risk_category <-  ifelse(score >= 5, "high_risk", 
                           ifelse(score >= 3.75, "medium_risk", "low_risk"))
  
  ### Descriptives ### 
  sum_low_risk[i] <- sum(risk_category == "low_risk")
  sum_medium_risk[i] <- sum(risk_category == "medium_risk")
  sum_high_risk[i] <- sum(risk_category == "high_risk")
  sum_outcome_observed[i] <- sum(Outcome_observed == 1)
  sum_outcome_true[i] <- sum(Outcome_true == 1)
  sum_high_risk_observed[i] <- sum(Outcome_observed[risk_category == "high_risk"] == 1)
  sum_high_risk_true[i] <- sum(Outcome_true[risk_category == "high_risk"] == 1)
  sum_medium_risk_observed[i] <- sum(Outcome_observed[risk_category == "medium_risk"] == 1)
  sum_medium_risk_true[i] <- sum(Outcome_true[risk_category == "medium_risk"] == 1)
  sum_low_risk_observed[i] <- sum(Outcome_observed[risk_category == "low_risk"] == 1)
  sum_low_risk_true[i] <- sum(Outcome_true[risk_category == "low_risk"] == 1)
  mean_age_low_risk_true = mean(age[risk_category == "low_risk" & Outcome_true == 1])
  mean_age_medium_risk_true = mean(age[risk_category == "medium_risk" & Outcome_true == 1])
  mean_age_high_risk_true = mean(age[risk_category == "high_risk" & Outcome_true == 1])
  
  ### Observed outcome ROC analysis ###
  # Calculate predicted probabilities 
  observed.beta <- glm(Outcome_observed ~ X1 + X2 + X3 + X4 + X5 + X6, family="binomial")$coef
  observed.pred <- logit.pred(observed.beta, cbind(X1, X2, X3, X4, X5, X6))
  
  # Observed Outcome ROC
  ROC_observed[i] <- mis_ROC(Outcome_observed, observed.pred, 0, 0)$auc
  
  # Observed Outcome - misclassification corrected ROC
  ROC_observed_corrected[i] <- mis_ROC(Outcome_observed, observed.pred, g0, g1)$auc
  
}

#### Results ####

# Save table 
ROC_results <- data.frame(
  ROC_observed_outcome = ROC_observed,
  ROC_observed_corrected = ROC_observed_corrected,
  sum_low_risk = sum_low_risk,
  sum_medium_risk = sum_medium_risk,
  sum_high_risk = sum_high_risk,
  sum_outcome_observed = sum_outcome_observed,
  sum_outcome_true = sum_outcome_true, 
  sum_high_risk_observed = sum_high_risk_observed, 
  sum_high_risk_true = sum_high_risk_true,
  sum_medium_risk_observed = sum_medium_risk_observed,
  sum_medium_risk_true = sum_medium_risk_true,
  sum_low_risk_observed = sum_low_risk_observed,
  sum_low_risk_true = sum_low_risk_true,
  mean_age_low_risk_true = mean_age_low_risk_true,
  mean_age_medium_risk_true = mean_age_medium_risk_true,
  mean_age_high_risk_true = mean_age_high_risk_true
)

write.csv(
  ROC_results, 
  file = paste0(tabfolder, "/Non-dif_low-med-high.csv"),
  row.names = F
)

cat("Simulation complete! Results saved in the output_tables folder as Non-dif_low-med-high.csv\n")

#### END ####
