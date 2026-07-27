#Install relevant R packages
install.packages("lme4", type = "source")
install.packages("languageserver")
install.packages("lmPerm")
# Load the package
library(lme4)
library(lmPerm)
#import data 
#C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Graph_Theory\schaef_400\flat_fc_dmn_all_within

#read csv file in of all possible fc connections from DMN to brain
fc_dmn_all_brain <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/fc_dmn_rest_brain_all_subjects.csv") # nolint # nolint: line_length_linter.

#read csv unique fc edges within dmn
upper_tri_within_dmn <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/upper_tri_fc_edges_within_dmn.csv")
# Display the first few rows of the data frame
dim(upper_tri_within_dmn) #returns dimensions of loaded data

subject_data <- read.csv("C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/subject_data.csv")
head(subject_data)

#need to select variables
overall_peak_drs_calc2 <- subject_data$overall_peak_drs_calc2 
nsqip_d <- subject_data$new_nsqip_d

#https://rdrr.io/github/cvoeten/permutes/man/perm.lmer.html

model_nsqipd <- lm(y ~ FC_edges + nsqip_d)
summary(model_nsqipd)

model_fcedges <- lm(y ~ upper_tri_within_dmn$FC_edge1 + nsqip_d)
summary(model_fcedges)

#attempt to join together the two dataframes
max_rows <- max(nrow(subject_data), nrow(upper_tri_within_dmn))


# Concatenate the extended data frames side-by-side
combined_df <- cbind(subject_data$overall_peak_drs_calc2, df2_extended)
combined_df_nsqipd <- cbind(overall_peak_drs_calc2, df2_extended, nsqip_d)

combined_df_nsqipd <- cbind(overall_peak_drs_calc2, fc_dmn_all_brain, nsqip_d)
## Attempt lm for loop and save results
# Create a list to store the results
results <- list()
results_perm <- list()
# Get the names of the predictor variables (excluding the response variable 'y')
predictor_vars <- names(combined_df)[names(combined_df) != "subject_data$overall_peak_drs_calc2"]

# Loop through each predictor variable
for (var in predictor_vars) {
  # Fit the linear model
  #model <- lm(as.formula(paste("subject_data$overall_peak_drs_calc2 ~", var)), data = combined_df)
   # nolint
  formula <- paste("subject_data$overall_peak_drs_calc2 ~", var)
  model_perm <- lmp(formula, data = combined_df, perm = "Exact") 

  # Save the summary of the model in the results list
  #results[[var]] <- summary(model)
  results_perm[[var]] <- summary(model_perm)
}

# Print the results
results

# permuted pvalues

pvalue <- results_perm$FC_edge3$p.values
print(pvalue)

pvalues_permute <- list()
for (var in predictor_vars) {
  pvalues_permute[[var]] <- results_perm[[var]]$p.values
}

#assessing specific results

coefficients <- results$FC_edge1$coefficients
print(coefficients)

predictor_vars <- names(combined_df_nsqipd)[names(combined_df_nsqipd) != "overall_peak_drs_calc2"]

# simple linear model predicting overall peak drs from individual fc edges within DMN
for (var in predictor_vars) {
## run general linear model
formula <- as.formula(paste("overall_peak_drs_calc2 ~", var, "+ nsqip_d"))
model <- lm(formula, data = combined_df_nsqipd) 
  # Save the summary of the model in the results list
  results[[var]] <- summary(model)
}

## Permute linear model
predictor_vars <- names(combined_df_nsqipd)[names(combined_df_nsqipd) != "overall_peak_drs_calc2"]
results_perm <- list()
p_values <- list()
# simple linear model predicting overall peak drs from individual fc edges within DMN
for (var in predictor_vars) {
  ## run general linear model
  formula <- as.formula(paste("overall_peak_drs_calc2 ~", var, "+ nsqip_d"))
  model_perm <- lmp(formula, data = combined_df_nsqipd, perm = "Exact", nperm = 1000) 
  # Save the summary of the model in the results list
  results_perm[[var]] <- summary(model_perm)
  p_values[[var]] <- model_perm$perm$P[2,] #saves permuted p-values for each fc edge in model
}

 #Convert the p_values list to a data frame
p_values_df <- as.data.frame(do.call(rbind, p_values))

# Count the number of p-values less than 0.05
num_significant <- sum(p_values_df$V1 < 0.05)
# estimate of coef edges
coef_edges <- list()
for (var in predictor_vars)
{
  coef_edges[[var]] <- results_perm[[var]]$coefficients[2,1]
  #print(coef_edges)
}
coef_edges_df <- as.data.frame(do.call(rbind, coef_edges))

#Find the coef_edges that are significant
# Add a column to coef_edges_df with the p-values
coef_edges_df$p_value <- p_values_df$V1

# Create a new data frame with only the rows where p_value < 0.05
significant_coef_edges_df <- coef_edges_df[coef_edges_df$p_value < 0.05, ]

# Remove the p_value column
significant_coef_edges_df$p_value <- NULL

# summary of negative beta coeffs
sum_negative_beta <- sum(significant_coef_edges_df$V1 < 0)

## save outputs from lm permuted model into csv for import to matlab
# Add a column to coef_edges_df with the p-values
coef_edges_df$p_value <- p_values_df$V1

# Create a new data frame with only the rows where p_value < 0.05
significant_coef_edges_df <- coef_edges_df[coef_edges_df$p_value < 0.05, ]

# Remove the p_value column

write.csv(significant_coef_edges_df, file = "C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/significant_coef_edges_df.csv", row.names = FALSE)

# Save significant_coef_edges_df to a CSV file
write.csv(, file = "significant_coef_edges.csv", row.names = FALSE)

# Save outputs of lm to a CSV file
write.csv(coef_edges_df, file = "C:/Users/natas/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/lm_model/permute_lm_within_dmn_coeffs_pval.csv", row.names = FALSE)



## create plots
library(ggplot2)
ggplot(significant_coef_edges_df, aes(x = rownames(significant_coef_edges_df), y = V1)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Variable", y = "Coefficient", title = "Significant Coefficients")

significant_coef_edges_df$sign <- ifelse(significant_coef_edges_df$V1 > 0, "positive", "negative")
# Create a bar plot of the significant coefficients with different colors for positive and negative coefficients
# Create a bar plot of the significant coefficients with different colors for positive and negative coefficients
ggplot(significant_coef_edges_df, aes(x = rownames(significant_coef_edges_df), y = V1, fill = sign)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("positive" = "blue", "negative" = "pink")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Variable", y = "Coefficient", title = "Significant Coefficients")










