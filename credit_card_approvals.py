import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV

# Load data
file_path = 'cc_approvals.data'
cc_apps = pd.read_csv(file_path, header=None)

# Replace '?' with NaN
cc_apps.replace('?', np.nan, inplace=True)

# Split features and target
X = cc_apps.iloc[:, :-1]
y = cc_apps.iloc[:, -1]

# Encode target
y = LabelEncoder().fit_transform(y)

# Identify numerical and categorical columns
numerical_mask = X.apply(lambda x: pd.to_numeric(x, errors='coerce')).notna().all()
numerical_columns = X.columns[numerical_mask].tolist()
categorical_columns = X.columns[~numerical_mask].tolist()

# Create preprocessing pipeline
preprocessor = ColumnTransformer(
    transformers=[
        ('num', SimpleImputer(strategy='mean'), numerical_columns),
        ('cat', Pipeline([
            ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
            ('onehot', OneHotEncoder(sparse=False, handle_unknown='ignore'))
        ]), categorical_columns)
    ])

# Create full pipeline
model = Pipeline([
    ('preprocessor', preprocessor),
    ('classifier', RandomForestClassifier(random_state=42))
])

# Split data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Parameter grid
param_grid = {
    'classifier__n_estimators': [100, 200],
    'classifier__max_depth': [10, 20, None]
}

# Grid search
grid_search = GridSearchCV(model, param_grid, cv=5, scoring='accuracy')
grid_search.fit(X_train, y_train)

# Store and print results
best_score = grid_search.best_score_
print(f"Best parameters: {grid_search.best_params_}")
print(f"Best cross-validation score: {best_score:.3f}")
print(f"Test set score: {grid_search.score(X_test, y_test):.3f}")