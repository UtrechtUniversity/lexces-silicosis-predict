#### About this script ####
## Title: Reconstructed diagnostic performance parameters for ILO 1/0 cutoff
##
## Author: Javier Mancilla Galindo
## ORCiD: https://orcid.org/0000-0002-0718-467X
##
## Purpose: This script calculates diagnostic performance measures for the 
## chest X-ray ILO 1/0 threshold against HRCT, using data reported in a 
## systematic review with meta-analysis:
## https://thorax.bmj.com/lookup/doi/10.1136/thorax-2024-BTSabstracts.23. 
## A plot is generated to model the trade-off between sensitivity and specificity, 
## using the calculated DOR.

## Note: This script is sourced into the main quarto markdown file (.qmd).
## Refer to README file of repository for instructions on how to use. 

#### START ####

# Load data from Howlett, et al. (2024) study available through: 
# https://github.com/pjhowlett/da_silic_cxr
HRCT <- readxl::read_excel(path = paste0(inputfolder,"/HRCT_1_1.xlsx"))

# Compute Sensitivity and Specificity for each study
HRCT$Sn <- HRCT$TP / (HRCT$TP + HRCT$FN)  # Sensitivity
HRCT$Sp <- HRCT$TN / (HRCT$TN + HRCT$FP)  # Specificity

# Estimate pooled sensitivity and specificity using random-effects model
HRCT_Sensitivity <- metaprop(HRCT$TP, HRCT$TP + HRCT$FN, 
                             common=FALSE, random=TRUE, 
                             sm="PLOGIT", studlab=HRCT$ID)

HRCT_Specificity <- metaprop(HRCT$TN, HRCT$TN + HRCT$FP, 
                             common=FALSE, random=TRUE, 
                             sm="PLOGIT", studlab=HRCT$ID)

# Extract pooled estimates and confidence intervals
pooled_Sn <- plogis(HRCT_Sensitivity$TE.random)
pooled_ci_Sn <- plogis(c(HRCT_Sensitivity$lower.random, HRCT_Sensitivity$upper.random))
pooled_Sp <- plogis(HRCT_Specificity$TE.random)
pooled_ci_Sp <- plogis(c(HRCT_Specificity$lower.random, HRCT_Specificity$upper.random))

# Obtain pooled DOR #
# log-DOR for each study; add 1 to TP, FP, FN, and TN to avoid zeros
HRCT$log_DOR <- log(((HRCT$TP + 1) * (HRCT$TN + 1)) / ((HRCT$FP + 1) * (HRCT$FN + 1)))

# Standard error of log DOR
HRCT$SE_log_DOR <- sqrt((1 / (HRCT$TP + 1)) + (1 / (HRCT$FP + 1)) + (1 / (HRCT$FN + 1)) + (1 / (HRCT$TN + 1)))

# Meta-analysis of DOR using random-effects model
HRCT_DOR <- metagen(
  TE = HRCT$log_DOR,       
  seTE = HRCT$SE_log_DOR,  
  studlab = HRCT$ID,      
  common = FALSE,     
  random = TRUE,      
  sm = "OR",
  level.ma = 0.95,
  exclude = which(HRCT$ID == "Talini 1995")
)

# The Talini study was excluded only for the calculation of the pooled DOR
# due to its large influence not allowing trade-off to converge with 95%CI
# from the meta-analysis

pooled_log_DOR <- HRCT_DOR$TE.random  
pooled_DOR <- exp(pooled_log_DOR)   
pooled_DOR_ci <- exp(c(HRCT_DOR$lower.random, HRCT_DOR$upper.random))

# Correction factor applied to DOR to better match possible combinations 
# of sensitivity and specificity in the trade-off plot 
pooled_DOR_corrected <- pooled_DOR + 15

#### Plot for the tradeoff of sensitivity and specificity #### 

# Define a sequence of Sensitivity values from 0.01 to 0.99
Sn_values <- seq(0.01, 0.99, length.out = 100)

# Compute corresponding Specificity values using pooled DOR. 
# The calculate_SP function is defined in environment if sourcing this script.
Sp_values <- calculate_Sp(Sn_values, pooled_DOR)
Sp_corrected <- calculate_Sp(Sn_values, pooled_DOR_corrected)
Sp_lower_curve <- calculate_Sp(Sn_values, pooled_DOR_ci[1])
Sp_upper_curve <- calculate_Sp(Sn_values, pooled_DOR_ci[2])

# Create data frame for plotting
plot_data <- data.frame(Sn_values, Sp_values, Sp_lower_curve, Sp_upper_curve)

# Generate the plot
plot_tradeoff_meta <- ggplot(plot_data, aes(x = Sn_values, y = Sp_values)) +
  geom_rect(aes(xmin = pooled_ci_Sn[1], xmax = pooled_ci_Sn[2], 
                ymin = pooled_ci_Sp[1], ymax = pooled_ci_Sp[2], fill = "Overlap"),
            alpha = 0.2) +
  geom_line(aes(color = "uncorrected"), linewidth = 1.2) +
  geom_line(aes(y = Sp_corrected, color = "corrected"), linewidth = 1.2) +
  geom_ribbon(aes(ymin = Sp_lower_curve, ymax = Sp_upper_curve, fill = "Trade-off"), alpha = 0.2) +
  geom_point(data = HRCT, aes(x = Sn, y = Sp), color = "black", size = 2, alpha = 0.6) +  
  scale_fill_manual(values = c("Overlap" = "burlywood1", "Trade-off" = "gray")) +
  scale_color_manual(values = c("uncorrected" = "darkred", "corrected" = "green4")) +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Sensitivity and Specificity trade-off (Howlett, et al. 2025)",
       x = "Sensitivity",
       y = "Specificity",
       tag = "B",
       fill = "95% CI",  
       color = "DOR",
       caption = paste("Corrected DOR = ", round(pooled_DOR_corrected))
       ) +
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

# Combine plots and save
ggsave(
  filename = file.path(figfolder, "Howlett-Hoy_trade-off_Sn_Sp.png"),
  plot = arrangeGrob(plot_tradeoff_Hoy, plot_tradeoff_meta, nrow = 2),
  width = 8, height = 10, dpi = 600, bg = "white"
)

#### END ####