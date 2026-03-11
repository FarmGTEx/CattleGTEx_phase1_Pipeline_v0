#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier

# ============================================================
# 1️ Command-line arguments
# ============================================================
if len(sys.argv) != 4:
    print("Usage: python train_breed_model.py <N> <selector> <method>")
    sys.exit(1)

N = int(sys.argv[1])          # e.g. 200
selector = sys.argv[2]        # e.g. "AED"
method = sys.argv[3].lower()  # e.g. "xgboost"

print(f"Job start | SNPs: {N} | Selector: {selector} | Method: {method}")

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

# ============================================================
# 3️  Load genotype matrix
# ============================================================
matrix_path = f"{selector}_{N}_ped.raw"
matrix = pd.read_csv(matrix_path, sep=r"\s+", header=0)

# drop PLINK metadata columns (FID IID PAT MAT SEX PHENOTYPE)
matrix = matrix.drop(columns=matrix.columns[1:6])
matrix = matrix.dropna(axis=1, how="any")

# set FID as index
matrix.index = matrix.iloc[:, 0]
matrix = matrix.drop(matrix.columns[0], axis=1)

# attach breed info
breed_map = dict(zip(all_sample["Sample"], all_sample["Breed"]))
matrix["Breed"] = matrix.index.map(breed_map)

# remove samples without breed label
matrix = matrix.dropna(subset=["Breed"])

y = matrix["Breed"].astype(str)
X = matrix.drop(columns=["Breed"])

print(f"Loaded genotype matrix: {X.shape[0]} samples × {X.shape[1]} SNPs")

# ============================================================
# 4️  Train/test split
# ============================================================
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=202406
)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# ============================================================
# 5️  Define model & Train
# ============================================================
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=202406)

# --- model definition ---
if method == "svp_p":
    model = SVC(kernel="linear", probability=False)
elif method == "svp_r":
    model = SVC(kernel="rbf", probability=False)
elif method == "rf":
    model = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=202406)
elif method == "knn":
    model = KNeighborsClassifier(n_neighbors=3, n_jobs=-1)
elif method == "xgboost":
    model = XGBClassifier(
        objective="multi:softmax",
        num_class=len(np.unique(y_train)), 
        eval_metric="merror",
        n_estimators=200,
        learning_rate=0.1,
        n_jobs=-1,
        random_state=202406,
        verbosity=0
    )
else:
    print(f"Error:  method '{method}'")
    sys.exit(1)

# --- Train and Predict ---
if method == "xgboost":

    le = LabelEncoder()
    y_train_enc = le.fit_transform(y_train)
    y_test_enc = le.transform(y_test)
   
    model.num_class = len(np.unique(y_train_enc))

    cv_scores = cross_val_score(model, X_train_scaled, y_train_enc, cv=cv, scoring="accuracy", n_jobs=-1)
    cv_acc = np.mean(cv_scores)
    cv_sd = np.std(cv_scores, ddof=1)
    print(f"CV accuracy (XGBoost): {cv_acc:.4f}")

    model.fit(X_train_scaled, y_train_enc)
    y_pred_enc = model.predict(X_test_scaled)
    y_pred = le.inverse_transform(y_pred_enc)
    y_true = y_test.reset_index(drop=True)

else:
    cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=cv, scoring="accuracy", n_jobs=-1)
    cv_acc = np.mean(cv_scores)
    cv_sd = np.std(cv_scores, ddof=1)
    print(f"CV accuracy: {cv_acc:.4f}")

    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    y_true = y_test.reset_index(drop=True)

# ============================================================

# ============================================================
# 6️  Evaluate & save results
# ============================================================
y_true = y_test.reset_index(drop=True)
test_acc = accuracy_score(y_true, y_pred)
cm = confusion_matrix(y_true, y_pred, labels=np.unique(y_true))

prefix = f"{method}_{selector}_{N}"
np.savetxt(f"{prefix}_cv_overall.txt", [cv_acc], fmt="%.5f")
np.savetxt(f"{prefix}_cv_sd.txt", [cv_sd], fmt="%.5f")
np.savetxt(f"{prefix}_leave_overall.txt", [test_acc], fmt="%.5f")

pd.DataFrame(cm, index=np.unique(y_true), columns=np.unique(y_true)).to_csv(
    f"{prefix}_leave_matrix.txt", sep="\t", index=True
)

report = classification_report(y_true, y_pred, digits=4)
with open(f"{prefix}_report.txt", "w") as f:
    f.write(report)

print(f"  - CV accuracy: {cv_acc:.4f}")
print(f"  - Test accuracy: {test_acc:.4f}")