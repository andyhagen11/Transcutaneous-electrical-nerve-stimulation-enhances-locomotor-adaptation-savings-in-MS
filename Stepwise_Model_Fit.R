# This script fits a regression model and assesses model fit for linear fixed effects or mixed effect models (depending on data type)
# Optionally performs stepwise regression 

# Load necessary libraries
library(lmerTest)
rm(list=setdiff(ls(), "data"))

# Print warnings as they come up for troubleshooting purposes
options(warn=1)

# Fit the model - These models specifically test the impact of different predictors on step length asymmetry (SLA) savings
#For fNIRS channel specific hypotheses
  formula = "EarlyADSave_SLA ~ S14.D14_EAd	+	S16.D7_EAd" 

#For fNIRS all ROIs
  #formula = "EarlyADSave_SLA ~ PMd +	PMv +	M1 + S1 +	SPL +	IPL"
  
#For demographic variables
  #formula = "EarlyADSave_SLA ~ Group + Condition + Sex + Age + FastLeg + EOFirm + ECFirm + EOFoam + ECFoam + Composite + SymbolDigit + CombinedVT + Avg_StepLength_B + StepLengthAsym_Normalized_B + BMI + RPE + Activity"

model = lm(formula, data = data, subset = ParticipantID != "SBNT_MS_007")

#Make model stepwise if desired
model = step(model)

# Create ANOVA Table
anovaResult = anova(model)

# Get summary of the model for coefficients and statistics (this includes R-squared, coefficients, etc.)
modelSummary = summary(model)

# Extract R-squared and Adjusted R-squared
rSquared = modelSummary$r.squared
adjRSquared = modelSummary$adj.r.squared

# Extract AIC and BIC
aicValue = AIC(model)
bicValue = BIC(model)

# Export coefficients (β coefficients)
coefficients = modelSummary$coefficients

# Extract confidence intervals for effect sizes
confIntervals = confint(model)

# Create output text with ANOVA results, model fit, and effect sizes
text = capture.output({
  cat("Model summary:\n")
  print(summary(model))  
  
  cat("ANOVA Table:\n")
  print(anovaResult)
  
  cat("\nModel Fit Statistics:\n")
  cat("R-squared: ", rSquared, "\n")
  cat("Adjusted R-squared: ", adjRSquared, "\n")
  cat("AIC: ", aicValue, "\n")
  cat("BIC: ", bicValue, "\n")
  
  cat("\nModel Coefficients (Effect Sizes):\n")
  print(coefficients)
  
  cat("\nConfidence Intervals for Coefficients:\n")
  print(confIntervals)
})

# Export the output to a text file
writeLines(text, file.path("/Users/andyhagen/Library/CloudStorage/OneDrive-Colostate/Colorado State/SNL/SBNT/Statistics/Correlations/Baseline Adj/SavingsfromROIs_Stepwise/SavingsfromROIs_Stepwise_MS_TENSOFF.txt"))

print(text)
