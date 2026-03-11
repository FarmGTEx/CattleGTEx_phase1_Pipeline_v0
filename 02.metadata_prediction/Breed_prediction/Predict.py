#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier

# =======================
# 1️ Command-line arguments
# =======================
#N = 10000
#selector ="PLSR"
#method = "svp_p"
N = int(sys.argv[1])         # Number of SNPs
selector = str(sys.argv[2])  # Feature selection method
method = str(sys.argv[3])    # Model type
print(f"Process {N} SNPs selected by {selector}, the method is {method}")

# ============================================================
# 2️  Load sample information
# ============================================================
rna_path = " " #metadata for RNA-seq samples
wgs_path = " " #metadata for WGS samples from reference panel

rna_df = pd.read_csv(rna_path, sep="\t", header=0, usecols=[0, 10])
wgs_df = pd.read_csv(wgs_path, sep="\t", header=0, usecols=[0, 2])
rna_df.columns = wgs_df.columns = ["Sample", "Breed"]

# naming normalization
rna_df["Sample"] = np.where(
    rna_df["Sample"].str.contains("lung|lymph|rumen", case=False),
    rna_df["Sample"] + ".bqsr",
    rna_df["Sample"] + "_bqsr"
)
all_sample = pd.concat([rna_df, wgs_df], ignore_index=True)

# =======================
# 3️  Load genotype matrix
# =======================

matrix = pd.read_csv(f"{selector}_{N}_ped.raw", sep=r"\s+", header=0)

matrix.index = matrix['IID']

matrix = matrix.drop(columns=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], errors='ignore')

matrix = matrix.dropna(axis=1)

matrix['breed'] = all_sample.set_index('Sample').reindex(matrix.index)['Breed']


# =======================
# 4️ Split samples into known and unknown sets
# =======================
train_data = matrix[matrix['breed'] != "Unknown"].copy()
train_data['breed'] = train_data['breed'].astype(str)

predict_data = matrix[matrix['breed'] == "Unknown"].copy()

X_train = train_data.drop(columns=['breed']).values
y_train = train_data['breed'].values
X_predict = predict_data.drop(columns=['breed']).values

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_predict_scaled = scaler.transform(X_predict)

# =======================
# 5  Define model
# =======================
if method == "svp_p":
    model = SVC(kernel="linear", probability=False, random_state=202406)
elif method == "svp_r":
    model = SVC(kernel="rbf", probability=False, random_state=202406)
elif method == "rf":
    model = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=202406)
elif method == "knn":
    model = KNeighborsClassifier(n_neighbors=3, n_jobs=-1)
elif method == "xgboost":
    le = LabelEncoder()
    y_train_enc = le.fit_transform(y_train)
    model = XGBClassifier(
        objective="multi:softmax",
        num_class=len(np.unique(y_train_enc)),
        eval_metric="merror",
        n_estimators=200,
        learning_rate=0.1,
        n_jobs=-1,
        random_state=202406,
        verbosity=0
    )
else:
    print(f"Error: Unknown method '{method}'")
    sys.exit(1)

# =======================
# 6️  Train model and Predict
# =======================
if method == "xgboost":
    model.fit(X_train_scaled, y_train_enc)
    y_pred_enc = model.predict(X_predict_scaled)
    y_pred = le.inverse_transform(y_pred_enc)
else:
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_predict_scaled)

# =======================
# 7️  Save result
# =======================
reference = pd.read_csv(" ", sep="\t")
reference['Predict_breed'] = reference['Breed_Class']

reference['Sample_id'] = reference.iloc[:,0].apply(
    lambda x: f"{x}.bqsr" if any(sub in x for sub in ["lung","lymph","rumen"]) else f"{x}_bqsr"
)

predict_indices = reference['Sample_id'].isin(predict_data.index)
reference.loc[predict_indices, 'Predict_breed'] = y_pred

reference.to_csv(f"reference_predict_{N}{method}.txt", sep="\t", index=False)
print(f"Final prediction saved to reference_predict_{N}{method}.txt")
