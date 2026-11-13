# -----------------------------------------------------------------------------
# NMWI.py 
# Function:
#   - Fit the NMWI logistic regression model.
#
# Auxiliary code fragment executed from Jupyter notebooks using:
#     %run -i Codes/Auxiliary_files/NMWI.py
#
# This file is not intended to be run as a standalone Python script. It relies on variables and imports already defined in the notebook
# environment which are made available through the interactive namespace via the -i option.

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score


C = 0.9375 # Best regularization parameters for relative abundance

# Model parameters
final_params = dict(
    random_state=42,
    penalty="l1",
    solver="liblinear",
    class_weight="balanced",
    C=C,
    max_iter=2000
)

NMWI = LogisticRegression(**final_params)
NMWI.fit(X_new, y.values)

y_train_pred = NMWI.predict(X_new)
scores_all = NMWI.decision_function(X_new)
train_bal_acc = balanced_accuracy_score(y.values, y_train_pred)

# NMWI coefficients
coef = NMWI.coef_.flatten()
feature_names = X_new.columns
coefficients = pd.DataFrame(coef, index=feature_names, columns=["Coefficient"])
sorted_coefficients = coefficients.sort_values("Coefficient", ascending=False)
NMWI_coefs = sorted_coefficients[(sorted_coefficients['Coefficient']>0) | (sorted_coefficients['Coefficient']<0)]
# NMWI_coefs.reset_index().to_csv('NMWI_coefficients.csv')

num_pos  = (sorted_coefficients["Coefficient"] > 0).sum()
num_neg  = (sorted_coefficients["Coefficient"] < 0).sum()
num_zero = (sorted_coefficients["Coefficient"] == 0).sum()
num_coef = len(sorted_coefficients)

NMWI_preds = X_new.reset_index()[['BioProject Number', 'Sample ID']]
NMWI_preds['NMWI'] = scores_all.tolist()
NMWI_preds.merge(sra_data[['Sample ID', 'Category', 'Disease']], on = 'Sample ID', how = 'left')

final_coef = sorted_coefficients[(sorted_coefficients['Coefficient']>0) | (sorted_coefficients['Coefficient']<0)]
print(f'Positive coefs: {final_coef.gt(0).sum().item()}')
print(f'Negative coefs: {final_coef.lt(0).sum().item()}')
print(f'zero coefs: {sorted_coefficients.eq(0).sum().item()}')
