# This code is saved into a script. 

# --------------------------------------------------
# PER-FOLD LOSO PERFORMANCE USING BEST C
# --------------------------------------------------

fold_bas = []
loso_loop = []
coef_rows = []
selected_sets = []
loso = LeaveOneGroupOut()
loso_scores = np.full(len(y), np.nan)

bio_ids = np.array(X_new.index.get_level_values(0))
y_arr = y.values.ravel()

bio_ids_arr = np.asarray(bio_ids)
fold_sizes = []
for train_idx, test_idx in loso.split(X_new, y_arr, groups=bio_ids):
    X_tr = X_new.iloc[train_idx]
    X_te = X_new.iloc[test_idx]
    y_tr = y_arr[train_idx]
    y_te = y_arr[test_idx]

    model = LogisticRegression(**final_params)
    model.fit(X_tr, y_tr)
    loso_scores[test_idx] = model.decision_function(X_te)

    y_hat = model.predict(X_te)
    ba = balanced_accuracy_score(y_te, y_hat)
    fold_bas.append(ba)
    fold_sizes.append(len(test_idx))

    proj_te = pd.Index(bio_ids_arr[test_idx]).unique().tolist()[0]

    loso_loop.append({
        "BioProject": proj_te,
        "n_test": len(test_idx),
        "pos_test": int(y_te.sum()),
        "neg_test": int((y_te == 0).sum()),
        "BA": ba
    })

    coefs = model.coef_.ravel()
    fold_name = pd.Index(bio_ids_arr[test_idx]).unique().tolist()[0]
    row = pd.Series(coefs, index=feature_names, name=fold_name)
    coef_rows.append(row)

    selected_sets.append(set(np.array(feature_names)[coefs != 0]))
    

coef_df = pd.DataFrame(coef_rows)
loso_df = pd.DataFrame(loso_loop).sort_values("BA")
mean_cv_ba = float(np.mean(fold_bas)) if fold_bas else np.nan
