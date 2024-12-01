import pandas as pd
from datetime import datetime
import os

# Define directory paths
base_dir = "GWAS_Project"
raw_data_path = os.path.join(base_dir, "data", "raw", "Harvard_Assessment_Chap2.csv")
cleaned_data_path = os.path.join(base_dir, "data", "cleaned", "Harvard_Assessment_Chap2_Cleaned.csv")
log_file_path = os.path.join(base_dir, "logs", "cleaning_log.txt")

# Ensure directories exist
os.makedirs(os.path.dirname(raw_data_path), exist_ok=True)
os.makedirs(os.path.dirname(cleaned_data_path), exist_ok=True)
os.makedirs(os.path.dirname(log_file_path), exist_ok=True)

# Load the dataset
df = pd.read_csv(raw_data_path)

# Create a unique patient ID
df['patient_id'] = range(1, len(df) + 1)

# Standardize Treatment column
df['treatment'] = df['Treatment A'].fillna('Treatment B').apply(lambda x: 0 if x == 'Treatment A' else 1)

# Standardize Age column
df['age'] = pd.to_numeric(df['Age of patient'], errors='coerce').replace({pd.NA: "Unknown"})

# Standardize Gender column
gender_map = {'Male': 0, 'Female': 1, 'M': 0, 'F': 1, '?': "Unknown"}
df['gender'] = df['Patient gender'].map(gender_map).fillna("Unknown")

# Standardize Height column
def standardize_height(height):
    if isinstance(height, str):
        if "'" in height and '"' in height:
            feet, inches = height.split("'")
            inches = inches.replace('"', '')
            return int(feet) * 12 + int(inches)
        elif '"' in height:
            return int(height.replace('"', ''))
        elif 'cm' in height.lower():
            return round(float(height.replace('cm', '')) / 2.54, 2)
    return "Unknown"

df['height_inches'] = df['Height'].apply(standardize_height)

# Standardize Tumor Stage column
tumor_map = {'I': 1, 'II': 2, 'III': 3, 'IV': 4}
df['tumor_stage'] = df['Tumor stage'].map(tumor_map).fillna("Unknown")

# Standardize Date Enrolled column
def parse_date(date):
    try:
        return (datetime.strptime(date, '%m/%d/%y') - datetime(1999, 1, 1)).days // 30
    except:
        return "Unknown"

df['months_since_start'] = df['Date enrolled'].apply(parse_date)

# Standardize Complications column
complications_map = {'no': 0, 'No': 0, 'yes': 1, 'Yes': 1, 'Y': 1, 'N': 0, '?': "Unknown"}
df['complications'] = df['Complications?'].map(complications_map).fillna("Unknown")

# Drop original columns for clarity
df = df[['patient_id', 'treatment', 'age', 'gender', 'height_inches', 'tumor_stage', 'months_since_start', 'complications']]

# Save cleaned dataset
df.to_csv(cleaned_data_path, index=False)

# Save a cleaning log
with open(log_file_path, "w") as log_file:
    log_file.write("Data cleaning completed successfully.\n")
    log_file.write(f"Raw data file: {raw_data_path}\n")
    log_file.write(f"Cleaned data file: {cleaned_data_path}\n")
    log_file.write(f"Columns cleaned: {', '.join(df.columns)}\n")

print(f"Data cleaned and saved to {cleaned_data_path}")
print(f"Log file saved to {log_file_path}")
