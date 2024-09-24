#### About this script ####
## Title: Estimation of the impact of outcome misclassification in a
##        pneumoconiosis diagnostic prediction rule
## Author: Javier Mancilla Galindo
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

# Number of repetitions
num_repetitions <- 5000

# Sample size for every new simulation
N <- 1291 # Chosen to be similar to diagnostic rule development sample size

# Frequency distribution of predictors in original OEM paper
prop_1 <- 690/1291 # Age > 40
prop_2 <- 642/1291 # Current smoker
prop_3 <- 868/1291 # High exposed job title
prop_4 <- 830/1291 # Work in the construction industry > 15 years
prop_5 <- 145/1291 # Feeling unhealthy
prop_6 <- 183/1291 # Standardized residual FEV1 <-1.0

# Prediction rule model effect size parameters in OEM paper
B0 <- -6.33
B1 <- 0.72 # Age > 40
B2 <- 0.70 # Current smoker
B3 <- 1.14 # High exposed job title
B4 <- 1.00 # Work in the construction industry > 15 years
B5 <- 0.846 # Feeling unhealthy
B6 <- 0.916 # Standardized residual FEV1 <-1.0

# Misclassification parameters were calculated from data reported in a 
# recent study in artificial stone benchtop industry workers
# https://onlinelibrary.wiley.com/doi/10.1111/resp.14755
Sn <- 0.575   # Sensitivity
Sp <- 0.971  # Specificity
g0 <- 1 - Sp  # gamma0 (false-positive rate)
g1 <- 1 - Sn  # gamma1 (false-negative rate)

# Initialize vectors to store results
ROC_observed <- vector(length = num_repetitions)
ROC_observed_corrected <- vector(length = num_repetitions)
sum_low_risk <- vector(length = num_repetitions)
sum_high_risk <- vector(length = num_repetitions)
sum_outcome_observed <- vector(length = num_repetitions)
sum_outcome_true <- vector(length = num_repetitions)
sum_high_risk_observed <- vector(length = num_repetitions)
sum_high_risk_true <- vector(length = num_repetitions)
sum_low_risk_observed <- vector(length = num_repetitions)
sum_low_risk_true <- vector(length = num_repetitions)

#### Simulation of data, scores, and AUC estimates for every new sample ####

for (i in 1:num_repetitions) {
  # Simulate Data
  X1 <- rbinom(N, 1, prop_1) # Age > 40
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
  Outcome_true <- misclassify_reverse(Outcome_observed, Sn, Sp)
  
  # Calculate scores according to the diagnostic prediction rule
  score <- X1*1 + X2*1 + X3*1.5 + X4*1.5 + X5*1.25 + X6*1.25
  
  # Categorize into high and low risk. 
  risk_category <- ifelse(score < 5, "low_risk", "high_risk")
  
  ### Descriptives ### 
  sum_low_risk[i] <- sum(risk_category == "low_risk")
  sum_high_risk[i] <- sum(risk_category == "high_risk")
  sum_outcome_observed[i] <- sum(Outcome_observed==1)
  sum_outcome_true[i] <- sum(Outcome_true==1)
  sum_high_risk_observed[i] <- sum(Outcome_observed[risk_category == "high_risk"] == 1)
  sum_high_risk_true[i] <- sum(Outcome_true[risk_category == "high_risk"] == 1)
  sum_low_risk_observed[i] <- sum(Outcome_observed[risk_category == "low_risk"] == 1)
  sum_low_risk_true[i] <- sum(Outcome_true[risk_category == "low_risk"] == 1)
  
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
  sum_high_risk = sum_high_risk,
  sum_outcome_observed = sum_outcome_observed,
  sum_outcome_true = sum_outcome_true, 
  sum_high_risk_observed = sum_high_risk_observed, 
  sum_high_risk_true = sum_high_risk_true, 
  sum_low_risk_observed = sum_low_risk_observed,
  sum_low_risk_true = sum_low_risk_true
)

### The following lines of code only work when sourced from main qmd file: 
write.csv(
  ROC_results, 
  file = paste0(tabfolder, "/ROC_results_scenario_1.csv"),
  row.names = F
  )

#### END ####
