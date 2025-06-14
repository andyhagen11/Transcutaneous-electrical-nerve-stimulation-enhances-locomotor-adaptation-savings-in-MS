# This script fits a linear mixed-effects mode, performs an anova test, and calculates estimated marginal means pairwise comparisons
# This model formula was used for all analyses, including step length asymmetry and fNIRS outcomes

# Load required packages
library(lmerTest)
library(emmeans)

rm(list=setdiff(ls(), "data"))

# Print warnings as they come up for troubleshooting purposes
options(warn=1)

# Prior to running, import data (as "data") in long format with each factor or variable as a separate column

# Fit the Linear Mixed-Effects Model (with variable of choice)
formula = "EarlyAD_SLA ~ Group * Condition * Visit + (1 | Participant)"
model = lmer(formula, data = data, subset = Participant != "SBNT_MS_007")

# Create ANOVA Table
anovaResult = anova(model)
#print(anovaResult)

# Calculate Estimated Marginal Means of Condition while averaging across Visit for each Group
emmCondition = emmeans(model, ~ Condition | Group)

# Calculate Estimated Marginal Means for Condition within each Group and Visit
emmConditionPairs = emmeans(model, ~ Condition | Group * Visit)

# Calculate Estimated Marginal Mean for Visit within each group (averages across Condition)
emmVisit = emmeans(model, ~ Visit | Group)

#Calculate Estimated Marginal Means for Group (both for Condition, and separated by Visit and Condition)
emmGroup = emmeans(model, ~ Group| Condition)
emmGroupPairs = emmeans(model, ~ Group| Condition * Visit)

# Pairwise Comparisons (note: multiple comparison corrections (using FDR) were done separately outside of this script)
pairwiseCondition = pairs(emmCondition)
pairwiseConditionPairs = pairs(emmConditionPairs)
pairwiseVisit = pairs(emmVisit)
pairwiseGroup = pairs(emmGroup)
pairwiseGroupPairs = pairs(emmGroupPairs)

# Calculate Cohen's D
cohensDCondition = as.data.frame(eff_size(emmCondition, sigma = sigma(model), edf = df.residual(model)))
cohensDConditionPairs = as.data.frame(eff_size(emmConditionPairs, sigma = sigma(model), edf = df.residual(model)))
cohensDVisit = as.data.frame(eff_size(emmVisit, sigma = sigma(model), edf = df.residual(model)))
cohensDGroup = as.data.frame(eff_size(emmGroup, sigma = sigma(model), edf = df.residual(model)))
cohensDGroupPairs = as.data.frame(eff_size(emmGroupPairs, sigma = sigma(model), edf = df.residual(model)))

# Store summary output
summaryCondition = summary(pairwiseCondition)
summaryCondition$Cohens_d = cohensDCondition$effect.size
summaryConditionPairs = summary(pairwiseConditionPairs)
summaryConditionPairs$Cohens_d = cohensDConditionPairs$effect.size
summaryVisit = summary(pairwiseVisit)
summaryVisit$Cohens_d = cohensDVisit$effect.size
summaryGroup = summary(pairwiseGroup)
summaryGroup$Cohens_d = cohensDGroup$effect.size
summaryGroupPairs = summary(pairwiseGroupPairs)
summaryGroupPairs$Cohens_d = cohensDGroupPairs$effect.size

# Capture output as text
anovaText = capture.output(print(anovaResult))
conditionPairsText = capture.output(print(summaryConditionPairs))
visitText = capture.output(print(summaryVisit))
conditionText = capture.output(print(summaryCondition))
groupText = capture.output(print(summaryGroup))
groupPairsText = capture.output(print(summaryGroupPairs))
combinedText = c( anovaText,
                "\nEstimated Margnial Means for Pairwise Comparisons by Group:\n", 
                groupPairsText, 
                "\nEstimated Marginal Means for Pairwise Comparisons by Condition:\n", 
                conditionPairsText, 
                "\nEstimated Marginal Means for Group:\n", 
                groupText, 
                "\nEstimated Marginal Means for Visit:\n", 
                visitText, 
                "\nEstimated Marginal Means for Condition:\n", 
                conditionText)

# Write to txt file
#writeLines(combinedText, file.path("/Users/andyhagen/Library/CloudStorage/OneDrive-Colostate/Colorado State/SNL/SBNT/Statistics/Savings/SLA/BaselineAdj/LME_EarlyAdapt.txt"))

# Print to console
cat(combinedText, sep = "\n")


