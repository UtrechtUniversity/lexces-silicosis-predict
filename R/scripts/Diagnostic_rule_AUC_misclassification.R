#### About this script ####
## Title: Estimation of the impact of outcome misclassification in a
##        pneumoconiosis diagnostic prediction rule
## Author: Javier Mancilla Galindo
## ORCiD: https://orcid.org/0000-0002-0718-467X
##
## Purpose: This script follows the procedure described by Zawitowski, et al. 
## https://onlinelibrary.wiley.com/doi/10.1002/sim.7260 to estimate the impact 
## of outcome misclassification in a pneumoconiosis diagnostic prediction rule
## https://oem.bmj.com/lookup/doi/10.1136/oem.2006.027904. However, the original
## paper estimates misclassification based on the TRUE outcome. Since the 
## diagnostic rule paper used a misclassified outcome (chest X-ray), I want to
## obtain the reverse misclassification (what is the probability of silicosis)
## if HRCT (reference test) had been used instead of CXR. For that, the code
## in Zawitowski's paper was adapted to simulate data compatible with the 
## summary sample data and prediction model equation in the OEM paper.
## A total of 5000 iterations are made to obtain stable mean AUC estimates. 
##
## Note: Functions by Zawitowski, et al. are found elsewhere. This script is
## sourced into the main quarto markdown file (.qmd) after function loading. 
## Refer to README file of repository for instructions on how to use. 

#### START ####

#### Set parameters ####

# Number of repetitions
num_repetitions <- 5000

# Sample size for every new matrix
N <- 1300 # Chosen to be similar to diagnostic rule development sample size

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

# Misclassification parameters from recent study in artificial stone benchtop
# industry workers https://onlinelibrary.wiley.com/doi/10.1111/resp.14755
g0 <- 0.0286  #gamma0 (false-positive rate)
g1 <- 0.425  # gamma1 (false-negative rate)

# Reverse misclassification parameters 
## Reverse misclassification parameters, calculated from data from a 
## recent study in artificial stone benchtop industry workers. Calculation 
## of these is found in the main Silicosis_diagnostic_rule.qmd file. 
rev_g0 <- 0.85185  # P(HRCT+ | CXR+), probability of HRCT+ given CXR+
rev_g1 <- 0.81928  # P(HRCT- | CXR-), probability of HRCT- given CXR-

# Initialize vectors to store results
ROC_observed <- vector(length = num_repetitions)
ROC_observed_corrected <- vector(length = num_repetitions)
ROC_reverse_misclassified <- vector(length = num_repetitions)

#### Simulation of data and AUC estimates for every new sample ####

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
  
  # Simulate observed outcome based on prediction rule probability
  Outcome_observed <- rbinom(N, 1, p)
  
  # Simulate reverse-misclassified outcomes 
  Y <- misclassify_reverse(Outcome_observed, rev_g0, rev_g1)
  
  # Observed Outcome 
  true.beta <- glm(Outcome_observed ~ X1 + X2 + X3 + X4 + X5 + X6, family="binomial")$coef
  true.pred <- logit.pred(true.beta, cbind(X1, X2, X3, X4, X5, X6))
  ROC_observed[i] <- mis_ROC(Outcome_observed, true.pred, 0, 0)$auc
  
  # Observed Outcome - misclassification corrected 
  ROC_observed_corrected[i] <- mis_ROC(Outcome_observed, true.pred, g0, g1)$auc
  
  #cor.beta = Misclassify_logistic_IWLS(
   # cbind(X1, X2, X3, X4, X5, X6), Outcome_observed, g0, g1)$coefficients
  #cor.pred = logit.pred(cor.beta, cbind(X1, X2, X3, X4, X5, X6))
  #ROC_observed_corrected[i] <- mis_ROC(Outcome_observed, cor.pred, g0, g1)$auc
  
  # Reverse-misclassified Outcome
  mis.rev.beta <- glm(Y ~ X1 + X2 + X3 + X4 + X5 + X6, family="binomial")$coef
  mis.rev.pred <- logit.pred(mis.rev.beta, cbind(X1, X2, X3, X4, X5, X6))
  ROC_reverse_misclassified[i] <- mis_ROC(Y, mis.rev.pred, 0, 0)$auc
}

#### Results ####

# Save table 
ROC_results <- data.frame(
  ROC_observed_outcome = ROC_observed,
  ROC_observed_corrected = ROC_observed_corrected,
  ROC_reverse_misclassified = ROC_reverse_misclassified
)

write.csv(
  ROC_results, 
  file = paste0(tabfolder, "/ROC_results_ILO1-1_Suarthana2007__ILO1-1_Hoy2024_reverse.csv"),
  row.names = F
  )

#### END ####
