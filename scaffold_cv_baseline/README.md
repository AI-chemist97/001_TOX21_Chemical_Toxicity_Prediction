# Tox21 12-target baseline

Chemical toxicity classification on all 12 Tox21 targets (MoleculeNet), evaluated with ROC-AUC on a scaffold split.

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
- Class imbalance (about 3-15% positives per target) handled with `class_weight="balanced_subsample"` (RF) or `scale_pos_weight` (LightGBM).
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

## Checking chemistry-based hypotheses

A model score does not show what the model relies on, so I picked four substructures with a chemical reason and checked whether their active rate is really higher in Tox21. I know quinones and hydroquinone as redox-active molecules because my master's research detected hydroquinone electrochemically. The hypotheses were written down before I looked at the numbers.

| Substructure | Hypothesis | Active rate with the group | Without |
|---|---|---|---|
| Quinone | SR-ARE (reactive oxygen species switch on the oxidative stress response) | 29 of 36 (81%) | 16% |
| Hydroquinone, catechol | SR-ARE (same reason) | 23 of 45 (51%) | 16% |
| Organotin | SR-MMP (damage to the mitochondrial membrane) | 11 of 16 (69%) | 16% |
| 3+ fused aromatic rings | NR-AhR (flat molecules bind the AhR receptor) | 55 of 82 (67%) | 11% |
| Nitroaromatic | SR-ARE (electrophilic stress) | 75 of 233 (32%) | 16% |

Hydroquinone itself (`Oc1ccc(O)cc1`) is in the data as one row and is active in SR-ARE and NR-AhR. NR-AhR was not one of my hypotheses and I have not explained it.

Code: `domain_analysis.py`. Full table: `domain_enrichment.csv`.

### Caveats

- Organotin has only 16 compounds, so its rate is unstable.
- I compared 11 groups against 12 targets, so some high values are probably chance. The table above only lists the hypotheses I fixed in advance.
- "Active" means a signal in a cell assay. It does not prove toxicity in humans, and cell death can also produce a signal.
- This is an association, not proof that the structure is the cause.

### Adding the flags to the model

I turned the 11 groups above into 0/1 columns, added them to the existing 2,265 features (2,276 in total) and ran the same scaffold 5-fold evaluation again. Code: `domain_features_cv.py`. Result: `domain_features_cv.csv`.

The mean ROC-AUC over the 12 targets was almost the same, 0.7945 before and 0.7937 after. Six targets went up a little and six went down a little. For all 12 targets the change was smaller than the fold-to-fold standard deviation, so I cannot say the flags helped. SR-MMP improved in all 5 folds, but only by +0.005.

My guess is that the Morgan fingerprint and RDKit descriptors already carry similar structural information, so the flags added nothing new. I have not checked this.
