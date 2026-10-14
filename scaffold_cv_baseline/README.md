# Tox21 12-target baseline

Chemical toxicity classification on all 12 Tox21 targets (MoleculeNet), evaluated with ROC-AUC on a scaffold split.

The code and this README in this folder were written with an AI coding tool (Claude Code). The original notebook in the repository root is the author's own work.

## Run

```bash
pip install -r requirements.txt
python tox21_baseline.py                  # RF + Morgan + RDKit descriptors, scaffold split
python tox21_baseline.py --split random   # same model, random split (expect higher scores)
python tox21_baseline.py --model lgbm --cv 5
```

The first run downloads `tox21.csv` (MoleculeNet / DeepChem mirror) into this folder.

## Method notes

- 12 targets: NR-AR, NR-AR-LBD, NR-AhR, NR-Aromatase, NR-ER, NR-ER-LBD, NR-PPAR-gamma, SR-ARE, SR-ATAD5, SR-HSE, SR-MMP, SR-p53.
- Missing labels (compound not tested on a target) are masked per target, not filled with 0.
- Class imbalance (about 3-16% positives per target) handled with `class_weight="balanced_subsample"` (RF) or `scale_pos_weight` (LightGBM).
- Scaffold split (Bemis-Murcko): no ring-system overlap between train and test, so the score reflects generalization to new scaffolds.
- Metrics: ROC-AUC (primary), PR-AUC (secondary, more informative at low positive rates).

## Results

Random Forest (300 trees, `class_weight="balanced_subsample"`), Morgan fingerprint (radius 2, 2048 bits) + 217 RDKit 2D descriptors = 2,265 features, 7,823 of 7,831 molecules parsed.
Scaffold-grouped 5-fold CV (`python tox21_baseline.py --model rf --split scaffold --features both --cv 5`), mean +- std over folds:

| Target | ROC-AUC | PR-AUC |
|---|---|---|
| NR-AR | 0.711 +- 0.144 | 0.372 +- 0.310 |
| NR-AR-LBD | 0.794 +- 0.129 | 0.438 +- 0.302 |
| NR-AhR | 0.824 +- 0.062 | 0.477 +- 0.194 |
| NR-Aromatase | 0.786 +- 0.041 | 0.292 +- 0.114 |
| NR-ER | 0.703 +- 0.115 | 0.410 +- 0.219 |
| NR-ER-LBD | 0.818 +- 0.046 | 0.399 +- 0.225 |
| NR-PPAR-gamma | 0.828 +- 0.042 | 0.323 +- 0.060 |
| SR-ARE | 0.766 +- 0.038 | 0.449 +- 0.128 |
| SR-ATAD5 | 0.824 +- 0.066 | 0.235 +- 0.085 |
| SR-HSE | 0.763 +- 0.049 | 0.309 +- 0.073 |
| SR-MMP | 0.883 +- 0.034 | 0.612 +- 0.227 |
| SR-p53 | 0.835 +- 0.030 | 0.311 +- 0.105 |
| **Mean of 12 targets** | **0.795** | **0.386** |

Notes:

- Targets with few positives (NR-AR, NR-AR-LBD, NR-ER) have a large fold-to-fold spread; read their means with care.
- No hyperparameter tuning was done; this is a baseline. The `--split random` comparison has not been run yet, so no random-split number is claimed here.
