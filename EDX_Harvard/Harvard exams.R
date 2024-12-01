source("https://bioconductor.org/biocLite.R")
biocLite("curatedOvarianData")
library("curatedOvarianData")

source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"))
data("GSE17260_eset")
mat.gene = exprs(GSE17260_eset)
y = mat.gene[1, ]  # Expression in the first column (target gene)
x = mat.gene[2:11, ]  # Expression values from columns 2 to 11 (predictors)
mod <- lm(y ~ t(x))  # Use transpose of x to match the dimensions
model_summary <- summary(mod)
r_squared <- model_summary$r.squared
print(r_squared)

# Prepare the response variable (vital_status)
y = GSE17260_eset$vital_status
y = ifelse(y == "deceased", 1, 0)

# Select the predictors: tumorstage and grade
x = cbind(as.numeric(GSE17260_eset$tumorstage), GSE17260_eset$grade)

# Combine into a data frame and remove missing values
bb = data.frame(y, x)
colnames(bb) = c("y", "tumorstage", "grade")  # Set appropriate column names

# Fit the logistic regression model
mod = glm(y ~ tumorstage + grade, data = bb, family = binomial(link = "logit"))

# Get predicted probabilities
pred_probs = predict(mod, type = "response")

# Compute the Brier Score
brier_score = mean((pred_probs - bb$y)^2)

# Print the Brier Score
print(brier_score)


install.packages("pROC")

# Load necessary library
library(pROC)

# Prepare the data
mat.gene = exprs(GSE17260_eset)
subset.gene = t(mat.gene[1:50, ])  # Use first 50 genes as predictors
y = GSE17260_eset$vital_status
y = ifelse(y == "deceased", 1, 0)

# Combine response and predictors into a data frame, removing NA values
uu = data.frame(na.omit(cbind(y, subset.gene)))

# Fit the logistic regression model
mod = glm(y ~ ., data = uu, family = binomial(link = "logit"))

# Get predicted probabilities
pred_probs = predict(mod, type = "response")

# Compute AUC using pROC
roc_curve = roc(uu$y, pred_probs)
auc_value = auc(roc_curve)

# Print the AUC
print(auc_value)


# Load necessary libraries
library(boot)

# Prepare the data
y = GSE17260_eset$vital_status
y = ifelse(y == "deceased", 1, 0)
x = cbind(as.numeric(GSE17260_eset$tumorstage), GSE17260_eset$grade)
bb = data.frame(na.omit(cbind(y, x)))

# Define the Brier score function
brier_score <- function(data, indices) {
  # Bootstrap resample
  boot_sample <- data[indices, ]
  y_boot <- boot_sample$y
  x_boot <- boot_sample[, -1]
  
  # Fit a logistic regression model
  mod <- glm(y_boot ~ ., data = data.frame(y = y_boot, x_boot), family = binomial(link = "logit"))
  
  # Get predicted probabilities
  pred_probs <- predict(mod, type = "response")
  
  # Calculate Brier score: mean squared error between observed and predicted
  brier <- mean((y_boot - pred_probs)^2)
  return(brier)
}

# Apply bootstrap with 1000 resamples
set.seed(123)
results <- boot(data = bb, statistic = brier_score, R = 1000)

# Compute the standard deviation of the Brier scores
sd_brier <- sd(results$t)
sd_brier






#######################
install.packages(c("tidyverse", "readr", "ggplot2", "dplyr", "psych"))
install.packages("Hmisc")  # For correlation and summary statistics
#########################

# Load the necessary library
library(readr)

# Read the TSV file
setwd("~/Desktop")
data <- read_tsv("04652-0001-Data.tsv")

# Display the first few rows of the data
head(data)


###Table 1
library(tidyverse)
library(readr)
data <- read.csv("path/to/MIDUS_data.csv")

data <- data %>%
  mutate(optimism = as.numeric(optimism),
         total_cholesterol = as.numeric(total_cholesterol),
         hdl = as.numeric(hdl),
         ldl = as.numeric(ldl),
         triglycerides = as.numeric(triglycerides))
colnames(data)

data <- data %>%
  mutate(optimism = as.numeric(B1SOPTIM), # Optimism score
         total_cholesterol = as.numeric(B1SA12C), # Total cholesterol
         cholesterol_frequency = as.numeric(B1SA12CY)) # Frequency of cholesterol RX in 30 days

summary(data$optimism)
summary(data$total_cholesterol)
summary(data$cholesterol_frequency)

library(ggplot2)

# Create the plot
plot <- ggplot(data, aes(x = optimism, y = total_cholesterol)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Optimism vs. Total Cholesterol",
       x = "Optimism",
       y = "Total Cholesterol")

# Save the plot as a JPEG file
ggsave("Optimism_vs_Total_Cholesterol.jpeg", plot = plot, width = 8, height = 6, dpi = 300)


####Figure 1
library(ggplot2)
library(dplyr)

# Define tertiles based on optimism scores
data <- data %>%
  mutate(optimism_tertile = case_when(
    optimism >= 6 & optimism <= 22 ~ "Lowest Tertile (6–22)",
    optimism >= 23 & optimism <= 26 ~ "Middle Tertile (23–26)",
    optimism >= 27 & optimism <= 30 ~ "Highest Tertile (27–30)"
  ))

# Create the histogram with tertile classification
figure_1 <- ggplot(data, aes(x = optimism, fill = optimism_tertile)) +
  geom_histogram(binwidth = 1, color = "black", position = "stack") +
  scale_fill_manual(values = c("black", "gray", "white"), name = "Optimism Tertiles") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 2)) + # Set x-axis range and breaks
  labs(
    title = "Figure 1: Frequency Distribution of Optimism Scores",
    x = "Optimism Score",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Display the plot
print(figure_1)

# Save the plot as a JPG file
ggsave("figure_1_optimism_distribution.jpg", plot = figure_1, width = 8, height = 6, dpi = 300)

library(ggplot2)
library(dplyr)

# Define tertiles based on optimism scores
data <- data %>%
  mutate(optimism_tertile = case_when(
    optimism >= 6 & optimism <= 22 ~ "Lowest Tertile (6–22)",
    optimism >= 23 & optimism <= 26 ~ "Middle Tertile (23–26)",
    optimism >= 27 & optimism <= 30 ~ "Highest Tertile (27–30)"
  ))

# Create the histogram with tertile classification
figure_1 <- ggplot(data, aes(x = optimism, fill = optimism_tertile)) +
  geom_histogram(binwidth = 1, color = "black", position = "stack") +
  scale_fill_manual(values = c("black", "gray", "white"), name = "Optimism Tertiles") +
  scale_x_continuous(limits = c(6, 30), breaks = seq(6, 30, by = 2)) + # Set x-axis range and breaks
  labs(
    title = "Figure 1: Frequency Distribution of Optimism Scores",
    x = "Optimism Score",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Display the plot
print(figure_1)

# Save the plot as a JPG file
ggsave("figure_1_optimism_distribution.jpg", plot = figure_1, width = 8, height = 6, dpi = 300)


