# Python
"""Tox21 12-target baseline: scaffold split, ROC-AUC per target, optional scaffold CV.

Usage:
    python tox21_baseline.py                      # RF, Morgan+descriptors, scaffold split
    python tox21_baseline.py --split random       # compare with a random split
    python tox21_baseline.py --model lgbm
    python tox21_baseline.py --cv 5               # scaffold-grouped 5-fold CV (mean +- std)
    python tox21_baseline.py --csv my_tox21.csv   # use a local file instead of downloading

Labels contain NaN (not every compound was tested on every target); those entries are masked
per target, never filled with 0.
"""
import argparse
import gzip
import shutil
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).parent
DATA_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz"
TARGETS = ["NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD",
           "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53"]
SEED = 42


# ---------- data ----------
def load_data(csv_path):
    """Return a dataframe with smiles + 12 target columns (NaN = not tested)."""
    path = Path(csv_path) if csv_path else HERE / "tox21.csv"
    if not path.exists():
        gz = HERE / "tox21.csv.gz"
        print("downloading", DATA_URL)
        urllib.request.urlretrieve(DATA_URL, gz)
        with gzip.open(gz, "rb") as f_in, open(path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    df = pd.read_csv(path)
    missing = [c for c in TARGETS + ["smiles"] if c not in df.columns]
    if missing:
        raise ValueError(f"missing columns: {missing}")
    return df


def clean_molecules(df):
    """Parse SMILES and drop rows RDKit cannot read. Returns (df, mols)."""
    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    keep = [m is not None for m in mols]
    print(f"parsed {sum(keep)}/{len(df)} molecules")
    df = df[keep].reset_index(drop=True)
    return df, [m for m in mols if m is not None]


# ---------- features ----------
def morgan_fp(mols, radius=2, n_bits=2048):
    out = np.zeros((len(mols), n_bits), dtype=np.uint8)
    for i, m in enumerate(mols):
        fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=n_bits)
        out[i, list(fp.GetOnBits())] = 1
    return out


def rdkit_descriptors(mols):
    """All RDKit 2D descriptors; inf/NaN -> 0 and clipped so one outlier cannot dominate."""
    names = [n for n, _ in Descriptors._descList]
    rows = []
    for m in mols:
        row = []
        for _, fn in Descriptors._descList:
            try:
                row.append(fn(m))
            except Exception:
                row.append(np.nan)
        rows.append(row)
    x = np.array(rows, dtype=np.float64)
    x[~np.isfinite(x)] = 0.0
    return np.clip(x, -1e6, 1e6), names


def build_features(mols, kind):
    parts = []
    if kind in ("fp", "both"):
        parts.append(morgan_fp(mols).astype(np.float32))
    if kind in ("desc", "both"):
        parts.append(rdkit_descriptors(mols)[0].astype(np.float32))
    return np.hstack(parts)


# ---------- splits ----------
def scaffold_of(mol):
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def scaffold_split(mols, frac_train=0.8, frac_valid=0.1):
    """Group by Bemis-Murcko scaffold; big scaffold groups go to train first.

    The same scaffold never appears on both sides, so the test set measures
    generalization to new ring systems rather than near-duplicates of training molecules.
    """
    groups = defaultdict(list)
    for i, m in enumerate(mols):
        groups[scaffold_of(m)].append(i)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), g[0]))
    n = len(mols)
    train, valid, test = [], [], []
    for g in ordered:
        if len(train) + len(g) <= frac_train * n:
            train += g
        elif len(valid) + len(g) <= frac_valid * n:
            valid += g
        else:
            test += g
    return np.array(train), np.array(valid), np.array(test)


def random_split(n, frac_train=0.8, frac_valid=0.1):
    idx = np.random.RandomState(SEED).permutation(n)
    a, b = int(frac_train * n), int((frac_train + frac_valid) * n)
    return idx[:a], idx[a:b], idx[b:]


# ---------- model ----------
def make_model(name, pos_weight):
    """Fresh classifier per target. pos_weight = n_neg / n_pos in that target's train rows."""
    if name == "rf":
        return RandomForestClassifier(n_estimators=300, class_weight="balanced_subsample",
                                      min_samples_leaf=2, n_jobs=-1, random_state=SEED)
    if name == "lgbm":
        from lightgbm import LGBMClassifier
        return LGBMClassifier(n_estimators=300, learning_rate=0.05, num_leaves=31,
                              scale_pos_weight=pos_weight, subsample=0.8, colsample_bytree=0.5,
                              random_state=SEED, verbose=-1)
    raise ValueError(name)


def eval_target(model_name, x, y, tr, te):
    """Train on labeled train rows of one target, score on labeled test rows. NaN labels masked."""
    tr = tr[~np.isnan(y[tr])]
    te = te[~np.isnan(y[te])]
    ytr, yte = y[tr].astype(int), y[te].astype(int)
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return None  # ROC-AUC is undefined with one class
    pos_weight = (ytr == 0).sum() / max((ytr == 1).sum(), 1)
    model = make_model(model_name, pos_weight).fit(x[tr], ytr)
    p = model.predict_proba(x[te])[:, 1]
    return {"n_train": len(tr), "n_test": len(te), "pos_rate_test": yte.mean(),
            "roc_auc": roc_auc_score(yte, p), "pr_auc": average_precision_score(yte, p)}


# ---------- runners ----------
def run_holdout(args, x, y_all, mols):
    split = scaffold_split(mols) if args.split == "scaffold" else random_split(len(mols))
    tr, va, te = split
    print(f"split={args.split}  train/valid/test = {len(tr)}/{len(va)}/{len(te)}")
    # TODO: use `va` to tune hyperparameters (Optuna); report `te` only once at the end.
    rows = []
    for t in TARGETS:
        r = eval_target(args.model, x, y_all[t].values, tr, te)
        if r:
            rows.append({"target": t, **r})
            print(f"{t:14s} ROC-AUC {r['roc_auc']:.3f}  PR-AUC {r['pr_auc']:.3f}  "
                  f"(test n={r['n_test']}, pos {r['pos_rate_test']:.1%})")
    return pd.DataFrame(rows)


def run_cv(args, x, y_all, mols):
    """Scaffold-grouped K-fold: every fold's test scaffolds are unseen in its train part."""
    groups = pd.factorize([scaffold_of(m) for m in mols])[0]
    rows = []
    for fold, (tr, te) in enumerate(GroupKFold(n_splits=args.cv).split(x, groups=groups)):
        for t in TARGETS:
            r = eval_target(args.model, x, y_all[t].values, tr, te)
            if r:
                rows.append({"target": t, "fold": fold, **r})
        print(f"fold {fold} done")
    df = pd.DataFrame(rows)
    return df.groupby("target")[["roc_auc", "pr_auc"]].agg(["mean", "std"]).round(3).reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="local tox21 csv (downloaded if omitted)")
    ap.add_argument("--features", choices=["fp", "desc", "both"], default="both")
    ap.add_argument("--model", choices=["rf", "lgbm"], default="rf")
    ap.add_argument("--split", choices=["scaffold", "random"], default="scaffold")
    ap.add_argument("--cv", type=int, default=0, help="K for scaffold-grouped CV (0 = single holdout)")
    ap.add_argument("--limit", type=int, default=0, help="use only the first N rows (quick test)")
    args = ap.parse_args()

    df = load_data(args.csv)
    if args.limit:
        df = df.head(args.limit)
    df, mols = clean_molecules(df)
    x = build_features(mols, args.features)
    print("feature matrix", x.shape)

    res = run_cv(args, x, df, mols) if args.cv else run_holdout(args, x, df, mols)
    out = HERE / f"results_{args.model}_{args.features}_{args.split if not args.cv else f'cv{args.cv}'}.csv"
    res.to_csv(out, index=False)
    if not args.cv:
        print(f"\nmean ROC-AUC over {len(res)} targets: {res['roc_auc'].mean():.3f}")
    else:
        print(res.to_string(index=False))
    print("saved", out)


if __name__ == "__main__":
    main()
