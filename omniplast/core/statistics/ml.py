"""Machine learning: bosques, SVM, boosting, clustering, validacion, interpretabilidad."""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.ensemble import (GradientBoostingClassifier, GradientBoostingRegressor,
                              RandomForestClassifier, RandomForestRegressor)
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score, mean_absolute_error,
                             mean_squared_error, precision_recall_curve, r2_score, roc_auc_score,
                             roc_curve, silhouette_score)
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


def _prep_xy(df, features, target, scale=True):
    data = df[features + [target]].apply(pd.to_numeric, errors="coerce").dropna()
    X, y = data[features].values, data[target].values
    if scale:
        X = StandardScaler().fit_transform(X)
    return X, y, data


def random_forest_classifier(df, features, target, test_size=0.3, cv=5, seed=42):
    X, y, data = _prep_xy(df, features, target, scale=False)
    classes, counts = np.unique(y, return_counts=True)
    if len(X) < 15 or len(classes) < 2 or counts.min() < 3:
        return {"modelo": "Random Forest", "error": "n≥15 y ≥3 casos por clase"}
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=test_size,
                                              stratify=y, random_state=seed)
    rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=2, random_state=seed,
                                class_weight="balanced").fit(X_tr, y_tr)
    pred = rf.predict(X_te)
    proba = rf.predict_proba(X_te)[:, 1] if len(classes) == 2 else None
    cv_scores = cross_val_score(rf, X, y, cv=StratifiedKFold(min(cv, counts.min()),
                                                             shuffle=True, random_state=seed))
    imp = permutation_importance(rf, X_te, y_te, n_repeats=30, random_state=seed)
    importances = pd.Series(imp.importances_mean, index=features).sort_values(ascending=False)
    metrics = {"accuracy_test": float(accuracy_score(y_te, pred)),
               "f1_test": float(f1_score(y_te, pred, average="weighted")),
               "cv_accuracy_media": float(cv_scores.mean()),
               "cv_accuracy_sd": float(cv_scores.std()),
               "n_train": len(y_tr), "n_test": len(y_te)}
    if proba is not None:
        metrics["roc_auc_test"] = float(roc_auc_score(y_te, proba))
    return {"modelo": "Random Forest (clasificación)", "metrics": metrics,
            "importancias": importances, "clases": classes.tolist(),
            "matriz_confusion": pd.DataFrame(confusion_matrix(y_te, pred),
                                             index=classes, columns=classes),
            "interpretacion": f"Exactitud test={metrics['accuracy_test']:.2%}, "
                              f"CV={metrics['cv_accuracy_media']:.2%}±{metrics['cv_accuracy_sd']:.2%}"}


def random_forest_regressor(df, features, target, test_size=0.3, seed=42):
    X, y, _ = _prep_xy(df, features, target, scale=False)
    if len(X) < 15:
        return {"modelo": "Random Forest", "error": "n≥15"}
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=test_size, random_state=seed)
    rf = RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=seed).fit(X_tr, y_tr)
    pred = rf.predict(X_te)
    imp = permutation_importance(rf, X_te, y_te, n_repeats=30, random_state=seed)
    importances = pd.Series(imp.importances_mean, index=features).sort_values(ascending=False)
    metrics = {"r2_test": float(r2_score(y_te, pred)),
               "rmse_test": float(np.sqrt(mean_squared_error(y_te, pred))),
               "mae_test": float(mean_absolute_error(y_te, pred)),
               "r2_train": float(rf.score(X_tr, y_tr))}
    return {"modelo": "Random Forest (regresión)", "metrics": metrics,
            "importancias": importances,
            "interpretacion": f"R² test={metrics['r2_test']:.3f}, RMSE={metrics['rmse_test']:.3g}"}


def svm_classifier(df, features, target, seed=42):
    X, y, _ = _prep_xy(df, features, target, scale=True)
    classes, counts = np.unique(y, return_counts=True)
    if len(X) < 15 or counts.min() < 3:
        return {"modelo": "SVM", "error": "n≥15 y ≥3 por clase"}
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, stratify=y, random_state=seed)
    svm = SVC(kernel="rbf", C=1.0, gamma="scale", probability=True, random_state=seed).fit(X_tr, y_tr)
    pred = svm.predict(X_te)
    return {"modelo": "SVM (RBF)", "metrics": {"accuracy": float(accuracy_score(y_te, pred)),
                                               "f1": float(f1_score(y_te, pred, average="weighted"))},
            "matriz_confusion": pd.DataFrame(confusion_matrix(y_te, pred), index=classes, columns=classes),
            "interpretacion": f"Exactitud={accuracy_score(y_te, pred):.2%}"}


def gradient_boosting(df, features, target, kind="auto", seed=42):
    X, y, _ = _prep_xy(df, features, target, scale=False)
    if len(X) < 15:
        return {"modelo": "Gradient Boosting", "error": "n≥15"}
    categorical = len(np.unique(y)) < 10 and np.all(np.equal(np.mod(y, 1), 0))
    if kind == "auto":
        kind = "clasificacion" if categorical else "regresion"
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, random_state=seed)
    if kind == "clasificacion":
        gb = GradientBoostingClassifier(n_estimators=300, learning_rate=0.05,
                                        max_depth=3, random_state=seed).fit(X_tr, y_tr)
        pred = gb.predict(X_te)
        out = {"modelo": "Gradient Boosting (clasificación)",
               "metrics": {"accuracy": float(accuracy_score(y_te, pred))}}
    else:
        gb = GradientBoostingRegressor(n_estimators=300, learning_rate=0.05,
                                       max_depth=3, random_state=seed).fit(X_tr, y_tr)
        pred = gb.predict(X_te)
        out = {"modelo": "Gradient Boosting (regresión)",
               "metrics": {"r2": float(r2_score(y_te, pred)),
                           "rmse": float(np.sqrt(mean_squared_error(y_te, pred)))}}
    imp = pd.Series(gb.feature_importances_, index=features).sort_values(ascending=False)
    out["importancias"] = imp
    return out


def dbscan_outliers(df, features, eps=1.5, min_samples=5):
    X = df[features].apply(pd.to_numeric, errors="coerce").dropna()
    Xs = StandardScaler().fit_transform(X)
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(Xs)
    labels = db.labels_
    out = df.loc[X.index].copy()
    out["cluster_dbscan"] = labels
    n_noise = int((labels == -1).sum())
    return {"tabla": out, "n_ruido": n_noise, "n_clusters": int(len(set(labels) - {-1})),
            "interpretacion": f"{n_noise} observaciones marcadas como atípicas espaciales "
                              f"({n_noise / len(X) * 100:.1f}%)" if len(X) else "sin datos"}


def roc_data(df, features, target, positive_class=None, seed=42):
    X, y, _ = _prep_xy(df, features, target, scale=True)
    classes = np.unique(y)
    if len(classes) != 2 or len(X) < 10:
        return {"modelo": "ROC", "error": "Se requiere objetivo binario con n≥10"}
    if positive_class is None:
        positive_class = classes[1]
    y_bin = (y == positive_class).astype(int)
    X_tr, X_te, y_tr, y_te = train_test_split(X, y_bin, test_size=0.3, random_state=seed)
    lr = LogisticRegression(max_iter=2000).fit(X_tr, y_tr)
    proba = lr.predict_proba(X_te)[:, 1]
    fpr, tpr, _ = roc_curve(y_te, proba)
    prec, rec, _ = precision_recall_curve(y_te, proba)
    return {"modelo": "Regresión logística + ROC", "fpr": fpr, "tpr": tpr,
            "auc": float(roc_auc_score(y_te, proba)), "precision": prec, "recall": rec,
            "odds_ratios": pd.Series(np.exp(lr.coef_[0]), index=features).sort_values(ascending=False),
            "interpretacion": f"AUC={roc_auc_score(y_te, proba):.3f} "
                              f"({'discriminación aceptable' if roc_auc_score(y_te, proba) > .7 else 'discriminación limitada'})"}


def pca_anomaly(df, features, contamination=None):
    """Deteccion de anomalias por error de reconstruccion PCA (stacked autoencoder-lite)."""
    X = df[features].apply(pd.to_numeric, errors="coerce").dropna()
    Xs = StandardScaler().fit_transform(X)
    n_comp = max(1, min(len(features) - 1, len(Xs) // 5))
    pca = PCA(n_components=n_comp).fit(Xs)
    recon = pca.inverse_transform(pca.transform(Xs))
    err = np.sqrt(((Xs - recon) ** 2).mean(axis=1))
    out = df.loc[X.index].copy()
    out["error_reconstruccion"] = err
    thr = np.percentile(err, 95)
    out["anomalia"] = err > thr
    return {"tabla": out, "umbral_p95": float(thr),
            "interpretacion": f"{int(out['anomalia'].sum())} observaciones con patrología "
                              f"multivariada atípica (error de reconstrucción > P95)"}


def silhouette_scan(df, features, k_min=2, k_max=8, scale=True):
    from sklearn.cluster import KMeans
    X = df[features].apply(pd.to_numeric, errors="coerce").dropna()
    if scale:
        X = StandardScaler().fit_transform(X)
    rows = []
    for k in range(k_min, min(k_max, len(X) - 1) + 1):
        km = KMeans(n_clusters=k, n_init=25, random_state=42).fit(X)
        sil = silhouette_score(X, km.labels_) if len(np.unique(km.labels_)) > 1 else np.nan
        rows.append({"k": k, "silueta": float(sil), "inercia": float(km.inertia_)})
    return pd.DataFrame(rows)
