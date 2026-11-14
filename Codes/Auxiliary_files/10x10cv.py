from pathlib import Path
import os, runpy, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, recall_score, confusion_matrix

warnings.filterwarnings("ignore")

# ---------------- Project / data ----------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.chdir(PROJECT_ROOT)
print(f"Working dir: {Path.cwd()}")

data = runpy.run_path(PROJECT_ROOT / "Codes/Auxiliary_files/Generate_data_index.py")
X_new = data["X"].copy()
y_arr = np.asarray(data["y"].values.ravel())

# ---------------- Model ----------------
C = 0.9375
final_params = {
    "random_state": 42,
    "penalty": "l1",
    "solver": "liblinear",
    "class_weight": "balanced",
    "C": C,
    "max_iter": 2000
}

# ---------------- Repeated 10x10 CV ----------------
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=10, random_state=42)
cv_splits = list(cv.split(X_new, y_arr))

# ---------------- Cutoffs ----------------
interval = 0.1
max_cutoff = 1.5
cutoffs = np.round(np.arange(0, max_cutoff + interval / 2, interval), 10)

oof_rows = []
cutoff_rows = []

for split_number, (train_idx, test_idx) in enumerate(cv_splits):
    repeat = split_number // 10 + 1
    fold = split_number % 10 + 1

    X_tr, X_te = X_new.iloc[train_idx], X_new.iloc[test_idx]
    y_tr, y_te = y_arr[train_idx], y_arr[test_idx]

    model_10 = LogisticRegression(**final_params)
    model_10.fit(X_tr, y_tr)

    decision_scores = model_10.decision_function(X_te)
    predicted_labels = model_10.predict(X_te)

    # Save individual OOF predictions
    for i, row_pos in enumerate(test_idx):
        oof_rows.append({
            "row_position": row_pos,
            "sample_index": X_new.index[row_pos],
            "repeat": repeat,
            "fold": fold,
            "true_label": y_te[i],
            "predicted_label": predicted_labels[i],
            "decision_score": decision_scores[i]
        })

    # Cutoff analysis
    for cutoff in cutoffs:
        keep = np.abs(decision_scores) >= cutoff
        n_retained = int(keep.sum())

        if n_retained == 0:
            continue

        y_keep = y_te[keep]
        y_pred_keep = predicted_labels[keep]

        bal_acc = balanced_accuracy_score(y_keep, y_pred_keep)
        recall = recall_score(y_keep, y_pred_keep, zero_division=0)

        tn, fp, fn, tp = confusion_matrix(
            y_keep, y_pred_keep, labels=[0, 1]
        ).ravel()

        tnr = tn / (tn + fp) if (tn + fp) > 0 else np.nan

        cutoff_rows.append({
            "cutoff": cutoff,
            "repeat": repeat,
            "fold": fold,
            "n_retained": n_retained,
            "samples_retained": n_retained / len(test_idx),
            "bal_acc": bal_acc,
            "recall": recall,
            "tnr": tnr
        })

# ---------------- OOF predictions ----------------
oof_predictions = pd.DataFrame(oof_rows).sort_values(
    ["repeat", "fold", "row_position"]
).reset_index(drop=True)

oof_predictions.to_csv(
    "Data/Output/10x10_oof_predictions.csv",
    index=False
)

# ---------------- Cutoff fold results ----------------
df_cutoff_fold = pd.DataFrame(cutoff_rows).sort_values(
    ["cutoff", "repeat", "fold"]
).reset_index(drop=True)

df_cutoff_fold.to_csv(
    "Data/Output/10x10_cutoff_fold_results.csv",
    index=False
)

# ---------------- Cutoff summary ----------------
df_cutoff_summary = (
    df_cutoff_fold
    .groupby("cutoff", as_index=False)
    .agg(
        samples_retained=("samples_retained", "mean"),
        samples_retained_sd=("samples_retained", lambda x: np.nanstd(x, ddof=1)),
        bal_acc=("bal_acc", "mean"),
        bal_acc_sd=("bal_acc", lambda x: np.nanstd(x, ddof=1)),
        recall=("recall", "mean"),
        tnr=("tnr", "mean")
    )
)

df_cutoff_summary.to_csv(
    "Data/Output/10x10_cutoff_summary.csv",
    index=False
)
