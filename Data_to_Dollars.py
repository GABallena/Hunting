import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import joblib

# Load and validate data
insurance = pd.read_csv('insurance.csv')
print("Data shape:", insurance.shape)
print("\nMissing values:\n", insurance.isnull().sum())

# Clean data
insurance['charges'] = pd.to_numeric(
    insurance['charges'].str.replace('[$,]', '', regex=True),
    errors='coerce'
)
insurance = insurance.dropna()

# Split features and target
X = insurance.drop('charges', axis=1)
y = insurance['charges']
numeric_features = ['age', 'bmi', 'children']
categorical_features = ['sex', 'smoker', 'region']

# Create robust preprocessor
preprocessor = ColumnTransformer([
    ('num', StandardScaler(), numeric_features),
    ('cat', OneHotEncoder(drop='first', sparse=False, handle_unknown='ignore'), 
     categorical_features)
])

# Create optimized pipeline
model = Pipeline([
    ('preprocessor', preprocessor),
    ('regressor', RandomForestRegressor(random_state=42))
])

# Split data with stratification
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Define comprehensive parameter grid
param_grid = {
    'regressor__n_estimators': [100, 200, 300],
    'regressor__max_depth': [None, 10, 20],
    'regressor__min_samples_split': [2, 5, 10]
}

# Perform parallel grid search
grid_search = GridSearchCV(
    model, param_grid, cv=5, scoring='r2', n_jobs=-1, verbose=1
)
grid_search.fit(X_train, y_train)

# Evaluate best model
best_model = grid_search.best_estimator_
test_r2_score = best_model.score(X_test, y_test)
cv_scores = cross_val_score(best_model, X, y, cv=5, scoring='r2')

# Print comprehensive results
print(f'\nModel Performance:')
print(f'Best Parameters: {grid_search.best_params_}')
print(f'R² Score (Test): {test_r2_score:.4f}')
print(f'R² Score (CV): {cv_scores.mean():.4f} (+/- {cv_scores.std()*2:.4f})')

# Save optimized model
joblib.dump(best_model, 'insurance_model.pkl')
print("\nModel saved successfully")