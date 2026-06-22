# -----------------------------------------------------------------------------
# NMWI.py 
# Function:
#   - Fit the NMWI logistic regression model.
#   - Evaluate model performance using repeated 10x10 stratified cross-validation (balanced accuracy)
#
# Auxiliary code fragment executed from Jupyter notebooks using:
#     %run -i Codes/Auxiliary_files/NMWI.py
#
# This file is not intended to be run as a standalone Python script. It relies on variables and imports already defined in the notebook
# environment which are made available through the interactive namespace via the -i option.

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

# 10 x 10 CV
y_arr = y.values.ravel()
model_10x10 = LogisticRegression(**final_params)

cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=10, random_state=42) 

cv_splits = []
cv_splits = list(cv.split(X_new, y_arr))

scores = cross_val_score(
    estimator=model_10x10,
    X=X_new,
    y=y.values,
    cv=cv,
    scoring='balanced_accuracy',
    n_jobs=-1
)

print('NMWI fitted; 10x10cv')
