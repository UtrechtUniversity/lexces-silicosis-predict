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
g1 <- 0.575  # gamma1 (true-positive rate)

# Initialize vectors to store results
ROC_true_results <- vector(length = num_repetitions)
ROC_misclassified_results <- vector(length = num_repetitions)

#### Simulation of data and AUC estimates for every new sample ####

for (i in 1:num_repetitions) {
  # Simulate Data
  X1 <- rbinom(N, 1, prop_1) # Age > 40
  X2 <- rbinom(N, 1, prop_2) # Current smoker
  X3 <- rbinom(N, 1, prop_3) # High exposed job title
  X4 <- rbinom(N, 1, prop_4) # Work in the construction industry > 15 years
  X5 <- rbinom(N, 1, prop_5) # Feeling unhealthy
  X6 <- rbinom(N, 1, prop_6) # Standardized residual FEV1 <-1.0
  
  p <- exp(
    B0 + B1*X1 + B2*X2 + B3*X3 + B4*X4 + B5*X5 + B6*X6) / 
    (1 + exp(B0 + B1*X1 + B2*X2 + B3*X3 + B4*X4 + B5*X5 + B6*X6)
     )
  T <- rbinom(N, 1, p)
  Y <- misclassify(T, g0, g1)
  
  # True Outcome Analysis
  true.beta <- glm(T ~ X1 + X2 + X3 + X4 + X5 + X6, family="binomial")$coef
  true.pred <- logit.pred(true.beta, cbind(X1, X2, X3, X4, X5, X6))
  ROC_true_results[i] <- mis_ROC(T, true.pred, 0, 0)$auc
  
  # Misclassified Outcome
  mis.beta <- glm(Y ~ X1 + X2 + X3 + X4 + X5 + X6, family="binomial")$coef
  mis.pred <- logit.pred(mis.beta, cbind(X1, X2, X3, X4, X5, X6))
  ROC_misclassified_results[i] <- mis_ROC(Y, mis.pred, 0, 0)$auc
}

#### Results ####

# Save table 
ROC_results <- data.frame(
  ROC_true = ROC_true_results,
  ROC_misclassified = ROC_misclassified_results
)

write.csv(
  ROC_results, 
  file = paste0(tabfolder, "/ROC_results_ILO1-1_Suarthana2007__ILO1-1_Hoy2024.csv"),
  row.names = F
  )

#### END ####
