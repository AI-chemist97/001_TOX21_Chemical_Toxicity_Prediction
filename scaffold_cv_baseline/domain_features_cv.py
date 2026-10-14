# Python
"""Does adding chemistry-motivated substructure flags improve scaffold-CV ROC-AUC?

Same model, folds and seeds as tox21_baseline.py (RF, Morgan + RDKit descriptors, scaffold-grouped
5-fold CV). The only change is 11 extra binary columns from domain_analysis.py.

Usage:
    python domain_features_cv.py        # writes domain_features_cv.csv
"""
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

import tox21_baseline as tb
from domain_analysis import PATTERNS, max_fused_aromatic_rings


def domain_flags(mols):
    cols = list(PATTERNS) + ["polycyclic_aromatic_3plus"]
    out = np.zeros((len(mols), len(cols)), dtype=np.float32)
    for i, m in enumerate(mols):
        for j, k in enumerate(PATTERNS):
            out[i, j] = any(m.HasSubstructMatch(p) for p in PATTERNS[k] if p is not None)
        out[i, -1] = max_fused_aromatic_rings(m) >= 3
    return out, cols


def main():
    df = tb.load_data(None)
    df, mols = tb.clean_molecules(df)
    x_base = tb.build_features(mols, "both")
    flags, cols = domain_flags(mols)
    x_dom = np.hstack([x_base, flags])
    print("features:", x_base.shape[1], "->", x_dom.shape[1], "flags:", cols)

    groups = pd.factorize([tb.scaffold_of(m) for m in mols])[0]
    rows = []
    for fold, (tr, te) in enumerate(GroupKFold(n_splits=5).split(x_base, groups=groups)):
        for t in tb.TARGETS:
            y = df[t].values
            rb = tb.eval_target("rf", x_base, y, tr, te)
            rd = tb.eval_target("rf", x_dom, y, tr, te)
            if rb and rd:
                rows.append({"target": t, "fold": fold, "auc_base": rb["roc_auc"], "auc_domain": rd["roc_auc"],
                             "pr_base": rb["pr_auc"], "pr_domain": rd["pr_auc"]})
        print("fold", fold, "done", flush=True)
    res = pd.DataFrame(rows)
    res["delta_auc"] = res["auc_domain"] - res["auc_base"]
    summary = res.groupby("target").agg(
        auc_base=("auc_base", "mean"), auc_domain=("auc_domain", "mean"),
        delta_mean=("delta_auc", "mean"), delta_std=("delta_auc", "std"),
        folds_improved=("delta_auc", lambda s: int((s > 0).sum()))).round(4).reset_index()
    summary.to_csv(tb.HERE / "domain_features_cv.csv", index=False)
    print(summary.to_string(index=False))
    print("mean over targets: base %.4f  domain %.4f  delta %.4f" % (
        summary.auc_base.mean(), summary.auc_domain.mean(), summary.delta_mean.mean()))


if __name__ == "__main__":
    main()
