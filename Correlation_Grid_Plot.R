#This script calculates pearson correlations and p values for multiple variables, then plots the data in a grid format

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra) 
rm(list=setdiff(ls(), "data"))


# Print warnings as they come up for troubleshooting purposes
options(warn=1)

# Prior to running, import data (as "data") in long format with each factor or variable as a separate column
data <- subset(data, ParticipantID != "SBNT_MS_007")

# Column names that need to be used, matching your dataset
variables = c("EarlyADSave_SLA",
               "S14.D14_EAd", "S16.D7_EAd")

# Convert necessary columns to numeric (if they aren't already)
data[variables] = lapply(data[variables], as.numeric)

# Filter the data by Group (MS and HC)
data_HC = data %>% filter(Group == "HC")
data_MS = data %>% filter(Group == "MS")

# Define a function to calculate correlations
correlation_test = function(df, variables) {
  results = data.frame(Variable = character(), Correlation = numeric(), P_Value = numeric())
  
  for (var in variables[-1]) {  # Exclude EarlyADSave_SLA from variable list
    cor_test = cor.test(df[[variables[1]]], df[[var]], method = "pearson", use = "complete.obs")
    results = rbind(results, data.frame(Variable = var, 
                                         Correlation = cor_test$estimate, 
                                         P_Value = cor_test$p.value))
  }
  return(results)
}

# Perform correlations for HC group
correlations_HC = correlation_test(data_HC, variables)
correlations_MS = correlation_test(data_MS, variables)

# Custom titles
hc_title = "HC Group fNIRS Correlations with EAdSave_SlA"
ms_title = "MS Group fNIRS Correlations with EAdSave_SlA"

# Convert data frames to character vectors (strings)
hc_content = capture.output(print(correlations_HC))
ms_content = capture.output(print(correlations_MS))

# Combine the title with the content and write
full_content = c(ms_title, ms_content,hc_title, hc_content)
#writeLines(full_content, file.path("/Users/andyhagen/Library/CloudStorage/OneDrive-Colostate/Colorado State/SNL/SBNT/Statistics/Correlations/SLA_Savings_TENSOFF.txt"))


# Function to create correlation plots with r values
create_plots = function(df, group_name) {
  plot_list = list()
  titles = c( "S14-D14 (SMA) at Early Adapt",
              "S16-D7 (SPL) at Early Adapt")
  pointColors = c("midnightblue",  "darkorange4")
  lineColors = c("blue", "darkorange")
  
  # Calculate correlations and store r values and p values
  correlations = correlation_test(df, variables)
  
  for (i in seq_along(variables[-1])) {  
    var = variables[-1][i] # Exclude EarlyADSave_SLA
    df_clean = df %>% select(all_of(variables)) %>% drop_na() # Remove rows with NA in relevant variables for plotting
    
    # Get the correlation value and p-value for the current variable
    r_value = round(correlations$Correlation[i], 3)
    p_value = round(correlations$P_Value[i], 4)
    
    plot = ggplot(df_clean, aes_string(x = var, y = variables[1])) +
      geom_point(size = 2, col = pointColors[i]) +
      geom_smooth(method = "lm", col = lineColors[i], se = TRUE) +  # Add regression line
      #ggtitle((titles[i])) +
      xlab("Beta") +
      ylab("SLA Savings") +
      scale_y_continuous(limits = c(-0.1, 0.1), breaks = seq(-0.08, 0.08, by = 0.04)) +  # Consistent scales
      theme_minimal(base_size = 14) 
      # + annotate("text", x = max(df_clean[[var]], na.rm = TRUE) * 0.5, 
      #          y = max(df_clean[[variables[1]]], na.rm = TRUE) * 1, 
      #          label = paste("r =", r_value, "\np =", p_value), 
      #          hjust = 0, vjust = 0, size = 4.25, color = lineColors[i])  # Adjusted hjust and vjust
    plot_list[[var]] = plot
  }
  return(plot_list)
}

# Create plots
plots_HC = create_plots(data_HC, "HC")
plots_MS = create_plots(data_MS, "MS")

# Arrange plots in a grid 
grid_HC = grid.arrange(grobs = plots_HC, ncol = 1)
grid_MS = grid.arrange(grobs = plots_MS, ncol = 1)
