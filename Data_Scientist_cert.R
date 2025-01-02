# First load your data
df <- read.csv("validation.csv")

# Count missing values in city column using multiple methods
missing_city <- sum(df$city == "--" | is.na(df$city) | df$city == "")

# Print result
missing_city





# Clean and standardize text data in the house_type column
df$house_type <- gsub("Det.", "Detached", df$house_type)
df$house_type <- gsub("Terr.", "Terraced", df$house_type)
df$house_type <- gsub("Semi", "Semi-detached", df$house_type)

# Clean area column by removing "sq.m." and converting to numeric
df$area <- as.numeric(gsub(" sq.m.", "", df$area))

# Standardize city names by trimming whitespace and converting to proper case
df$city <- tools::toTitleCase(trimws(df$city))

# Replace missing values
df$city[df$city == "" | df$city == "--" | is.na(df$city)] <- "Unknown"
df$house_type[df$house_type == "" | is.na(df$house_type)] <- "Unknown"
df$area[is.na(df$area)] <- mean(df$area, na.rm = TRUE)

# Print cleaned data summary
print("Unique house types:")
print(table(df$house_type))
print("\nUnique cities:")
print(table(df$city))




# Aggregate numeric variables: sale_price, area, months_listed by bedrooms
numeric_agg <- aggregate(
    cbind(sale_price, area, months_listed) ~ bedrooms, 
    data = clean_data,
    FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x))
)

# Aggregate categorical variables: house_type and city by bedrooms
cat_agg_house <- table(clean_data$bedrooms, clean_data$house_type)
cat_agg_city <- table(clean_data$bedrooms, clean_data$city)

# Aggregate dates: count sales by month
clean_data$month <- format(as.Date(clean_data$sale_date), "%Y-%m")
date_agg <- table(clean_data$month)

# Print results
print("Numeric Aggregations by Bedrooms:")
print(numeric_agg)
print("\nHouse Type Distribution by Bedrooms:")
print(cat_agg_house)
print("\nCity Distribution by Bedrooms:")
print(cat_agg_city)
print("\nSales by Month:")
print(date_agg)





# Install and load required packages
if (!require(caret)) install.packages("caret")
if (!require(randomForest)) install.packages("randomForest")
library(caret)
library(randomForest)

# Prepare data for modeling
model_data <- df[, !(names(df) %in% c('house_id', 'sale_date'))]
model_data$sale_price <- as.numeric(model_data$sale_price)

# Convert categorical variables to factors
categorical_cols <- sapply(model_data, is.character)
model_data[categorical_cols] <- lapply(model_data[categorical_cols], as.factor)

# Split data into training and testing sets
set.seed(42)
train_index <- createDataPartition(model_data$sale_price, p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Define cross-validation settings
ctrl <- trainControl(method = "cv", number = 5)

# Train multiple models
# 1. Random Forest
rf_model <- train(
    sale_price ~ .,
    data = train_data,
    method = "rf",
    trControl = ctrl,
    ntree = 100
)

# 2. Linear Regression
lm_model <- train(
    sale_price ~ .,
    data = train_data,
    method = "lm",
    trControl = ctrl
)

# Compare models
models <- list(RandomForest = rf_model, LinearRegression = lm_model)
results <- resamples(models)
summary(results)

# Use best model for final predictions
best_model <- rf_model  # Random Forest typically performs better
predictions <- predict(best_model, newdata = test_data)

# Calculate performance metrics
rmse <- sqrt(mean((test_data$sale_price - predictions)^2))
r2 <- cor(test_data$sale_price, predictions)^2

cat("Model Performance:\n")
cat("RMSE:", rmse, "\n")
cat("R-squared:", r2, "\n")
