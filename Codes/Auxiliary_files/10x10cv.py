from pathlib import Path
import os
import pandas as pd
import numpy as np
import glob, warnings, sys, runpy
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut, GridSearchCV, RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_curve, roc_auc_score

warnings.filterwarnings("ignore")


# ---------------- Project / data ----------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.chdir(PROJECT_ROOT)
print(f"Working dir: {Path.cwd()}")

data = runpy.run_path(PROJECT_ROOT / "Codes/Auxiliary_files/Generate_data_index.py")
X_new = data["X"].copy()
y_arr = np.asarray(data["y"].values.ravel())
C_grid = np.logspace(-2, 1, 40)


# ---------------- Model ----------------
model_params = dict(
    random_state=42,
    penalty="l1",
    solver="liblinear",
    class_weight="balanced",
    max_iter=2000
)

def select_C_10fold(X, y, seed=42):
    inner = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
    mean_ba = []

    for C in C_grid:
        fold_ba = []

        for tr, va in inner.split(X, y):
            model = LogisticRegression(**model_params, C=C)
            model.fit(X.iloc[tr], y[tr])
            fold_ba.append(
                balanced_accuracy_score(y[va], model.predict(X.iloc[va]))
            )

        mean_ba.append(np.mean(fold_ba))

    best_idx = np.argmax(mean_ba)
    return C_grid[best_idx], mean_ba[best_idx]


outer = RepeatedStratifiedKFold(
    n_splits=10,
    n_repeats=10,
    random_state=42
)

outer_results, pred_rows = [], []

for i, (tr, te) in enumerate(outer.split(X_new, y_arr)):
    repeat = i // 10 + 1
    fold = i % 10 + 1

    X_tr, X_te = X_new.iloc[tr], X_new.iloc[te]
    y_tr, y_te = y_arr[tr], y_arr[te]

    # Tune C using outer-training samples only
    best_C, inner_ba = select_C_10fold(X_tr, y_tr, seed=42+i)

    # Refit on complete outer-training set
    model = LogisticRegression(**model_params, C=best_C)
    model.fit(X_tr, y_tr)

    y_hat = model.predict(X_te)
    scores = model.decision_function(X_te)

    outer_results.append({
        "repeat": repeat,
        "fold": fold,
        "best_C": best_C,
        "inner_BA": inner_ba,
        "outer_BA": balanced_accuracy_score(y_te, y_hat),
        "outer_AUROC": roc_auc_score(y_te, scores)
    })

    # Save every held-out prediction
    for pos, idx, yt, yp, score in zip(te, X_new.index[te], y_te, y_hat, scores):
        sample_id = "|".join(map(str, idx)) if isinstance(idx, tuple) else str(idx)

        pred_rows.append({
            "row_position": pos,
            "sample_index": sample_id,
            "repeat": repeat,
            "fold": fold,
            "true_label": int(yt),
            "predicted_label": int(yp),
            "decision_score": score,
            "best_C": best_C
        })

    print(f"Repeat {repeat}, fold {fold} | C={best_C:.4g}")


# ---------------- Results ----------------
full_10x10_df = pd.DataFrame(outer_results)
oof_10x10 = pd.DataFrame(pred_rows)

# Mean balanced accuracy: average folds within repeat, then repeats
repeat_ba = full_10x10_df.groupby("repeat")["outer_BA"].mean()

# AUROC: pool all held-out predictions within each repeat
repeat_auc = oof_10x10.groupby("repeat").apply(
    lambda d: roc_auc_score(d["true_label"], d["decision_score"]),
    include_groups=False
)

repeat_summary = pd.DataFrame({
    "repeat": np.arange(1, 11),
    "balanced_accuracy": repeat_ba.values,
    "AUROC": repeat_auc.values
})

# print(f"\nfull 10x10 CV balanced accuracy: {repeat_ba.mean():.3f}")
# print(f"full 10x10 CV AUROC: {repeat_auc.mean():.3f}")


# ---------------- Export ----------------
OUTDIR = Path("Data")
OUTDIR.mkdir(exist_ok=True)

full_10x10_df.to_csv(
    OUTDIR / "10x10_fold_results.csv",
    index=False
)

oof_10x10.to_csv(
    OUTDIR / "10x10_oof_predictions.csv",
    index=False
)

repeat_summary.to_csv(
    OUTDIR / "10x10_repeat_summary.csv",
    index=False
)
