#### About this script ####
## Title: Reconstructed diagnostic performance parameters for ILO 1/1 cutoff
## against HRCT for the diagnosis of silicosis.
##
## Author: Javier Mancilla Galindo
## ORCiD: https://orcid.org/0000-0002-0718-467X
##
## Purpose: This script calculates diagnostic performance measures for the 
## chest X-ray ILO 1/1 threshold against HRCT, using data reported in a 
## systematic review with meta-analysis:
## https://thorax.bmj.com/lookup/doi/10.1136/thorax-2024-BTSabstracts.23. 
## The studies included in the meta-analysis are for an ILO 1/0 cutoff. Thus, 
## individual studies were re-examined and a total of 4 studies from which 
## diagnostic performance for the ILO 1/1 cutoff could be reconstructed 
## are included in this meta-analysis. 
## A plot is generated to model the trade-off between sensitivity and specificity, 
## using the pooled diagnostic odds ratio (DOR).

## Note: This script is sourced into the main quarto markdown file (.qmd).
## Refer to README file of repository for instructions on how to use. 

#### START ####

# Load data 
HRCT <- readxl::read_excel(path = paste0(inputfolder,"/HRCT_1_1.xlsx"))

# Compute Sensitivity and Specificity for each study
HRCT$Sn <- HRCT$TP / (HRCT$TP + HRCT$FN)  # Sensitivity
HRCT$Sp <- HRCT$TN / (HRCT$TN + HRCT$FP)  # Specificity

#### Sensitivity and Specificity Meta-analysis ####

# Estimate pooled sensitivity and specificity using random-effects model, 
# following Howlett, et al.: https://github.com/pjhowlett/da_silic_cxr 
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

#### Obtain pooled DOR ####

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
  level.ma = 0.95
)

pooled_log_DOR <- HRCT_DOR$TE.random  
pooled_DOR <- exp(pooled_log_DOR)   
pooled_DOR_ci <- exp(c(HRCT_DOR$lower.random, HRCT_DOR$upper.random))

#### Calculate Specificity from DOR #### 

# Define function
calculate_Sp <- function(Sn, DOR) {
  return(1 - (Sn / (Sn + ((1 - Sn) * DOR))))
}

# Define a sequence of Sensitivity values from 0.01 to 0.99
Sn_values <- seq(0.01, 0.99, length.out = 100)

# Compute corresponding Specificity values using pooled DOR. 
Sp_values <- calculate_Sp(Sn_values, pooled_DOR)
Sp_lower_curve <- calculate_Sp(Sn_values, pooled_DOR_ci[1])
Sp_upper_curve <- calculate_Sp(Sn_values, pooled_DOR_ci[2])

#### Plot #### 

# Create data frame for plotting
plot_data <- data.frame(Sn_values, Sp_values, Sp_lower_curve, Sp_upper_curve)

# Generate the plot
plot_tradeoff_meta <- ggplot(plot_data, aes(x = Sn_values, y = Sp_values)) +
  geom_rect(aes(xmin = pooled_ci_Sn[1], xmax = pooled_ci_Sn[2], 
                ymin = pooled_ci_Sp[1], ymax = pooled_ci_Sp[2], fill = "Overlap"),
            alpha = 0.2) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_ribbon(aes(ymin = Sp_lower_curve, ymax = Sp_upper_curve, fill = "Trade-off"), alpha = 0.2) +
  geom_point(data = HRCT, aes(x = Sn, y = Sp), color = "black", size = 2, alpha = 0.6) +  
  scale_fill_manual(values = c("Overlap" = "burlywood1", "Trade-off" = "gray")) +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Sensitivity and Specificity trade-off (CXR ILO 1/1 vs HRCT)",
       x = "Sensitivity",
       y = "Specificity",
       fill = "95% CI",  
       caption = paste("Diagnostic Odds Ratio (DOR) = ", round(pooled_DOR))
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
  filename = file.path(figfolder, "trade-off_Sn_Sp.png"),
  width = 8, height = 5, dpi = 600, bg = "white"
)

#### END ####