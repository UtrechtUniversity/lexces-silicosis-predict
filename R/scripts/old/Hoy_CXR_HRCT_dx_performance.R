#### About this script ####
## Title: Reconstructed diagnostic performance parameters for ILO 1/1 cutoff
##
## Author: Javier Mancilla Galindo
## ORCiD: https://orcid.org/0000-0002-0718-467X
##
## Purpose: This script calculates diagnostic performance measures for the 
## chest X-ray ILO 1/1 threshold against HRCT, using data reported in a 
## recent study in artificial stone benchtop industry workers
## https://onlinelibrary.wiley.com/doi/10.1111/resp.14755. A plot is 
## generated to model the trade-off between sensitivity and specificity, 
## using the calculated DOR.

## Note: This script is sourced into the main quarto markdown file (.qmd).
## Refer to README file of repository for instructions on how to use. 

#### START ####

# Input values from recreated 2x2 table 
TP <- 23 # True positives
TN <- 68 # True negatives 
FP <- 2  # False positives
FN <- 17 # False negatives 

# Calculate point estimate 
Sn <- TP/(TP+FN)
Sp <- TN/(FP+TN)
FPR <- 1-Sp
FNR <- 1-Sn

# Confidence intervals 
ci_Sn <- prop.test(TP, TP+FN, conf.level = 0.95)$conf.int
ci_Sp <- prop.test(TN, FP+TN, conf.level = 0.95)$conf.int
ci_FPR <- 1 - ci_Sp
ci_FNR <- 1 - ci_Sn

# Additional diagnostic performance measures 
PPV <- TP/(TP+FP)
NPV <- TN/(FN+TN)
LRpos <- Sn/(1-Sp)
LRneg <- (1-Sn)/Sp
Accuracy <- (TP+TN)/(TP+FP+FN+TN)

# Diagnostic Odds Ratio (DOR)
conf_matrix <- matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE)
# Compute the exact odds ratio and confidence intervals
results_DOR <- oddsratio.midp(conf_matrix, conf.level = 0.95)

DOR <- results_DOR$measure[2, 1]  # DOR estimate
DOR_lower <- results_DOR$measure[2, 2]  # Lower bound
DOR_upper <- results_DOR$measure[2, 3]  # Upper bound

#### Calculate Specificity from DOR #### 
calculate_Sp <- function(Sn, DOR) {
  return(1 - (Sn / (Sn + ((1 - Sn) * DOR))))
}

#### Plot for the tradeoff of sensitivity and specificity #### 

# Define a sequence of Sensitivity values from 0.01 to 0.99
Sn_values <- seq(0.01, 0.99, length.out = 100)

# Compute corresponding Specificity values using the DOR formula
Sp_values <- calculate_Sp(Sn_values, DOR)
Sp_values_lower <- calculate_Sp(Sn_values, DOR_lower)
Sp_values_upper <- calculate_Sp(Sn_values, DOR_upper)

# Create data frame for plotting
plot_data <- data.frame(Sn_values, Sp_values, Sp_values_lower, Sp_values_upper)

# Generate the plot

plot_tradeoff_Hoy <- ggplot(plot_data, aes(x = Sn_values, y = Sp_values)) +
  geom_rect(aes(xmin = ci_Sn[1], xmax = ci_Sn[2], ymin = ci_Sp[1], ymax = ci_Sp[2], fill = "Overlap"),
            alpha = 0.2) +
  geom_line(linewidth = 1.2, color = "cadetblue4") +
  geom_ribbon(aes(ymin = Sp_values_lower, ymax = Sp_values_upper, fill = "Trade-off"), alpha = 0.2) +
  ylim(0, 1) +
  scale_fill_manual(values = c("Overlap" = "burlywood1", "Trade-off" = "gray")) +
  labs(title = "Sensitivity and Specificity trade-off (based on Hoy, et al. 2024)",
       x = "Sensitivity",
       y = "Specificity",
       tag = "A",
       fill = "95% CI",  
       caption = "The Overlap region represents all possible combinations within 95% CIs, without assuming a trade-off.") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(colour = "black"),
    plot.margin = margin(5, 5, 5, 5, "mm"),
    plot.tag = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold")
  ) 

#### END ####