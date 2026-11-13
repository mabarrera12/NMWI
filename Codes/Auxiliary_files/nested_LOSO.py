from pathlib import Path
import os, runpy, warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

warnings.filterwarnings("ignore")


# ---------------- Project / data ----------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
os.chdir(PROJECT_ROOT)
print(f"Working dir: {Path.cwd()}")

data = runpy.run_path(PROJECT_ROOT / "Codes/Auxiliary_files/Generate_data_index.py")
X_new = data["X"].copy()
y_arr = np.asarray(data["y"].values.ravel())

groups = np.asarray(X_new.index.get_level_values(0))
feature_names = X_new.columns
classes = np.unique(y_arr)

if len(classes) != 2:
    raise ValueError("Binary outcome required.")

neg_label, pos_label = classes

OUTDIR = PROJECT_ROOT / "Data"
OUTDIR.mkdir(parents=True, exist_ok=True)


# ---------------- Model ----------------
model_params = dict(
    random_state=42,
    penalty="l1",
    solver="liblinear",
    class_weight="balanced",
    max_iter=2000
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

    sen, spe = div(tp, tp + fn), div(tn, tn + fp)
    pre = div(tp, tp + fp)
    f1 = div(2 * tp, 2 * tp + fp + fn)

    ba = balanced_accuracy_score(y_true, y_pred)
    strict_ba = np.nan if np.isnan(sen) or np.isnan(spe) else (sen + spe) / 2

    auc = np.nan
    if scores is not None and len(np.unique(y_true)) == 2:
        auc = roc_auc_score((y_true == pos_label).astype(int), scores)

    return {
        "Accuracy": ba,
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

            model = LogisticRegression(**model_params, C=C)
            model.fit(X.iloc[tr], y[tr])

            scores.append(
                balanced_accuracy_score(
                    y[va],
                    model.predict(X.iloc[va])
                )
            )

        mean_score = np.mean(scores) if scores else -np.inf

        if mean_score > best_score:
            best_C, best_score = C, mean_score

    return best_C, best_score


# ---------------- Nested outer LOSO ----------------
outer = LeaveOneGroupOut()
outer_splits = list(outer.split(X_new, y_arr, groups))

outer_pred = np.empty(len(y_arr), dtype=y_arr.dtype)
outer_scores = np.full(len(y_arr), np.nan)

rows, pred_rows, coef_rows = [], [], []

print(f"\nRunning {len(outer_splits)} outer LOSO folds")

for fold, (tr, te) in enumerate(outer_splits, 1):
    study = np.unique(groups[te])[0]

    Xtr, Xte = X_new.iloc[tr], X_new.iloc[te]
    ytr, yte = y_arr[tr], y_arr[te]

    # Select C using only outer-training studies
    best_C, inner_ba = select_C(Xtr, ytr, groups[tr])

    # Fit outer model
    model = LogisticRegression(**model_params, C=best_C)
    model.fit(Xtr, ytr)

    pred = model.predict(Xte)
    score = model.decision_function(Xte)

    outer_pred[te] = pred
    outer_scores[te] = score

    print(f"{fold}/{len(outer_splits)}: {study} | C={best_C:.4g}")

    # Fold-level metrics
    rows.append({
        "Row": "FOLD",
        "Fold": fold,
        "BioProject": study,
        "n_test": len(te),
        "best_C": best_C,
        "inner_LOSO_BA": inner_ba,
        **metrics(yte, pred, score)
    })

    # Per-sample OOF predictions
    for idx, yt, yp, sc in zip(X_new.index[te], yte, pred, score):
        pred_rows.append({
            "Sample_Index": idx,
            "BioProject": study,
            "Fold": fold,
            "True_label": yt,
            "Predicted_label": yp,
            "Prediction_score": sc
        })

    # Per-fold model coefficients
    coefs = model.coef_.ravel()
    coef_rows.append(
        pd.Series(coefs, index=feature_names, name=study)
    )


# ---------------- Results ----------------
results = pd.DataFrame(rows)
oof_df = pd.DataFrame(pred_rows)

coef_df = pd.DataFrame(coef_rows)
coef_df.index.name = "BioProject"

metric_cols = [
    "Accuracy",
    "Strict_BA",
    "AUROC",
    "Sensitivity",
    "Specificity",
    "Precision",
    "F1"
]

# Unweighted mean across held-out studies
mean_row = {
    "Row": "MEAN ACROSS FOLDS",
    "Fold": "",
    "BioProject": "",
    "n_test": results["n_test"].sum(),
    "best_C": np.nan,
    "inner_LOSO_BA": results["inner_LOSO_BA"].mean(),
    **{m: results[m].mean(skipna=True) for m in metric_cols}
}

# Pooled OOF performance across all samples
pooled_row = {
    "Row": "POOLED OUT-OF-FOLD",
    "Fold": "",
    "BioProject": "",
    "n_test": len(y_arr),
    "best_C": np.nan,
    "inner_LOSO_BA": np.nan,
    **metrics(y_arr, outer_pred, outer_scores)
}

final_df = pd.concat(
    [results, pd.DataFrame([mean_row, pooled_row])],
    ignore_index=True
)

# ---------------- Export ----------------
oof_df.to_csv(
    OUTDIR / "Output/nested_loso_oof_predictions.tsv",
    sep="\t", index=False, na_rep="NA"
)

final_df.to_csv(
    OUTDIR / "Output/nested_loso_summary.tsv",
    sep="\t", index=False, na_rep="NA"
)