#### About this script ####
## Title: Sample characteristics for simulations
## Author: Javier Mancilla Galindo
## ORCiD: https://orcid.org/0000-0002-0718-467X
##
## Purpose: This script sets the sample characteristics for the simulations
## by replicating construction workers' data from the original OEM paper.
##
## Note: This script is sourced into the main quarto markdown file (.qmd).
## Refer to README file of repository for instructions on how to use. 

# Sample size for every new simulation
N <- 1291 # Chosen to be similar to diagnostic rule development sample size

# Prediction rule model effect size parameters in OEM paper
B0 <- -6.33
B1 <- 0.72 # Age > 40
B2 <- 0.70 # Current smoker
B3 <- 1.14 # High exposed job title
B4 <- 1.00 # Work in the construction industry > 15 years
B5 <- 0.846 # Feeling unhealthy
B6 <- 0.916 # Standardized residual FEV1 <-1.0

# Frequency distribution of predictors in original OEM paper
prop_1 <- 690/1291 # Age > 40
prop_2 <- 642/1291 # Current smoker
prop_3 <- 868/1291 # High exposed job title
prop_4 <- 830/1291 # Work in the construction industry > 15 years
prop_5 <- 145/1291 # Feeling unhealthy
prop_6 <- 183/1291 # Standardized residual FEV1 <-1.0

# Age distribution
n1 <- 1254     # Outcome-negative participants
n2 <- 37       # Outcome-positive participants
mean1 <- 41.3  # Mean age for outcome-negative participants
mean2 <- 46.1  # Mean age for outcome-positive participants
sd1 <- 7.7
sd2 <- 7.9
# Weighted mean
age_mean <- ((n1 * mean1) + (n2 * mean2)) / (n1 + n2)
# Pooled standard deviation
age_sd <- sqrt((((n1 - 1) * sd1^2) + ((n2 - 1) * sd2^2)) / (n1 + n2 - 2))
