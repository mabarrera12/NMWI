# Execute from the terminal
from pathlib import Path
import os
import pandas as pd
import numpy as np
import glob, warnings, sys, runpy
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut, GridSearchCV, RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")
# Run only once to get the working directory
print("Working dir:", Path.cwd())
vars_dict = runpy.run_path("Codes/Auxiliary_files/Generate_data_index.py")
globals().update(vars_dict)
X_new = X.copy()

# Data/model setup
groups = np.asarray(X_new.index.get_level_values(0))
y_arr = np.asarray(y.values.ravel())
neg_label, pos_label = np.unique(y_arr)

model_params = dict(
    random_state=42, penalty="l1", solver="liblinear",
    class_weight="balanced", max_iter=2000
)

C_grid = np.logspace(-2, 1, 40)

def div(a, b):
    return np.nan if b == 0 else a / b

def metrics(y_true, y_pred, scores=None):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)

    tp = np.sum((y_true == pos_label) & (y_pred == pos_label))
    tn = np.sum((y_true == neg_label) & (y_pred == neg_label))
    fp = np.sum((y_true == neg_label) & (y_pred == pos_label))
    fn = np.sum((y_true == pos_label) & (y_pred == neg_label))

    sen = div(tp, tp + fn)
    spe = div(tn, tn + fp)
    pre = div(tp, tp + fp)
    f1 = div(2 * tp, 2 * tp + fp + fn)

    accuracy = balanced_accuracy_score(y_true, y_pred)
    strict_ba = np.nan if np.isnan(sen) or np.isnan(spe) else (sen + spe) / 2

    auc = np.nan
    if scores is not None and len(np.unique(y_true)) == 2:
        auc = roc_auc_score((y_true == pos_label).astype(int), scores)

    return {
        "Accuracy": accuracy,
        "Strict_BA": strict_ba,
        "AUROC": auc,
        "Sensitivity": sen,
        "Specificity": spe,
        "Precision": pre,
        "F1": f1
    }

def select_C(X, y, groups):
    logo = LeaveOneGroupOut()
    best_C, best_score = None, -np.inf

    for C in C_grid:
        scores = []
        for tr, va in logo.split(X, y, groups):
            if len(np.unique(y[tr])) < 2:
                continue
            m = LogisticRegression(**model_params, C=C)
            m.fit(X.iloc[tr], y[tr])
            scores.append(balanced_accuracy_score(y[va], m.predict(X.iloc[va])))

        mean_score = np.mean(scores) if scores else -np.inf
        if mean_score > best_score:
            best_C, best_score = C, mean_score

    return best_C, best_score

# Nested outer LOSO
outer = LeaveOneGroupOut()
outer_pred = np.empty(len(y_arr), dtype=y_arr.dtype)
outer_scores = np.full(len(y_arr), np.nan)
rows = []

for tr, te in outer.split(X_new, y_arr, groups):
    Xtr, Xte = X_new.iloc[tr], X_new.iloc[te]
    ytr, yte = y_arr[tr], y_arr[te]

    best_C, inner_ba = select_C(Xtr, ytr, groups[tr])

    model = LogisticRegression(**model_params, C=best_C)
    model.fit(Xtr, ytr)

    pred = model.predict(Xte)
    score = model.decision_function(Xte)

    outer_pred[te] = pred
    outer_scores[te] = score
    print(np.unique(groups[te])[0])

    rows.append({
        "Row": "FOLD",
        "BioProject": np.unique(groups[te])[0],
        "n_test": len(te),
        "best_C": best_C,
        **metrics(yte, pred, score)
    })

results = pd.DataFrame(rows)

metric_cols = [
    "Accuracy", "Strict_BA", "AUROC",
    "Sensitivity", "Specificity", "Precision", "F1"
]

# Mean across folds
mean_row = {
    "Row": "MEAN ACROSS FOLDS",
    "BioProject": "",
    "n_test": results["n_test"].sum(),
    "best_C": np.nan,
    **{m: results[m].mean(skipna=True) for m in metric_cols}
}

# Pooled out-of-fold
pooled_row = {
    "Row": "POOLED OUT-OF-FOLD",
    "BioProject": "",
    "n_test": len(y_arr),
    "best_C": np.nan,
    **metrics(y_arr, outer_pred, outer_scores)
}

final_df = pd.concat(
    [results, pd.DataFrame([mean_row, pooled_row])],
    ignore_index=True
)

print(final_df.round(3).to_string(index=False))