# Tox21 Toxicity Prediction

Rev. 0 | Created: 2026-10-12 | Updated: 2026-10-12 20:15 KST

Tox21 데이터의 12개 독성 타깃을 모두 예측하고, scaffold 기준으로 나눈 검증에서 ROC-AUC를 보고하는 baseline 프로젝트이다. 이전 notebook (NR-AR 한 타깃, 정확도 중심)을 12개 타깃과 ROC-AUC 중심의 평가로 다시 설계했다.

## 1. Purpose

### 1.1 Problem

화학 구조에서 독성 여부를 예측하는 분류 문제이다. 한 화합물이 모든 타깃에서 측정된 것은 아니고, 양성 비율이 낮아 정확도만으로는 성능을 판단할 수 없다.

### 1.2 Changes from the First Notebook

Table 1. First notebook versus this baseline

| Item | First notebook | This baseline |
|---|---|---|
| Targets | NR-AR 1개 | 12개 전체 |
| Imbalance | 다수 클래스 undersampling (118건 사용) | class weight, 데이터 전부 사용 |
| Split | 명시 없음 | Bemis-Murcko scaffold 기준 5-fold |
| Metric | 정확도 (50%, 67%) | ROC-AUC, PR-AUC |
| Missing labels | `dropna()` | 타깃별 masking |

첫 notebook의 문제는 약 24건짜리 test에서 정확도 67%라는 숫자가 우연과 구분되지 않았다는 점이다. 이 baseline은 그 평가 방식 자체를 바꾼 것이다.

## 2. Method

### 2.1 Data and Features

- 데이터: MoleculeNet Tox21, 7,831건 중 RDKit이 읽은 7,823건
- 특징: Morgan fingerprint (반경 2, 2,048 bit) + RDKit 2D descriptor 217개, 총 2,265개

### 2.2 Validation

- Scaffold 분할: 같은 Bemis-Murcko scaffold가 train과 test에 동시에 들어가지 않는다. 새로운 고리 구조에 대한 일반화를 잰다.
- 5-fold scaffold-grouped CV: fold마다 test scaffold는 train에 없다.
- 결측 label은 0으로 채우지 않고 타깃별로 제외했다.

### 2.3 Model

Random Forest (300 trees, `class_weight="balanced_subsample"`, `min_samples_leaf=2`, seed 42). 타깃마다 따로 학습했다.

## 3. Results

Table 2. Scaffold-grouped 5-fold CV, mean +- std over folds

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
| **Mean of 12** | **0.795** | **0.386** |

SR-MMP (0.883)와 SR-p53 (0.835) 등 양성이 많은 타깃에서 성능이 높고, 양성이 적은 NR-AR (0.711), NR-ER (0.703)은 낮고 fold 간 편차가 크다.

## 4. Limitations

+ Hyperparameter tuning을 하지 않은 baseline이다.
+ 양성이 적은 타깃의 평균은 fold에 따라 크게 달라진다. NR-AR의 ROC-AUC 표준편차는 0.144이다.
+ Random split과의 비교는 아직 실행하지 않았다. 그 비교 결과가 나오기 전에는 scaffold 분할이 얼마나 더 어려운지 수치로 말하지 않는다.
+ 그래프 신경망, 사전학습 분자 표현은 사용하지 않았다.

## 5. Next Steps

+ Random split과 scaffold split의 차이를 같은 모델로 측정한다.
+ LightGBM과 비교하고 validation 구간에서 hyperparameter를 조정한다.
+ 타깃 간 상관을 이용하는 multi-task 모델을 시험한다.

## 6. Run

```bash
pip install -r requirements.txt
python tox21_baseline.py --model rf --split scaffold --features both --cv 5
```

첫 실행에서 `tox21.csv`를 자동으로 내려받는다. 결과는 `results_rf_both_cv5.csv`에 저장된다.

## 화학 지식으로 확인해 본 것

모델 점수만 보면 모델이 무엇을 보고 판단하는지 알 수 없어서, 화학적으로 이유가 있는 구조를 먼저 정하고 Tox21 데이터에서 활성 비율이 정말 높은지 확인했습니다. 석사 연구에서 하이드로퀴논을 전기화학적으로 검출했기 때문에 퀴논과 하이드로퀴논이 산화·환원을 반복하는 물질이라는 점을 알고 있었고, 이를 바탕으로 가설을 세웠습니다. 가설은 숫자를 보기 전에 적어 두었습니다.

| 구조 | 가설 | 구조가 있는 물질의 활성 비율 | 없는 물질 |
|---|---|---|---|
| 퀴논 | SR-ARE 활성 (활성산소가 생겨 산화 스트레스 반응이 켜짐) | 36개 중 29개 (81%) | 16% |
| 하이드로퀴논, 카테콜 | SR-ARE 활성 (같은 이유) | 45개 중 23개 (51%) | 16% |
| 유기주석 | SR-MMP 활성 (미토콘드리아 막에 손상) | 16개 중 11개 (69%) | 16% |
| 3환 이상 방향족 | NR-AhR 활성 (평평한 분자가 AhR 수용체에 결합) | 82개 중 55개 (67%) | 11% |
| 나이트로 방향족 | SR-ARE 활성 (친전자성 스트레스) | 233개 중 75개 (32%) | 16% |

데이터에 하이드로퀴논(`Oc1ccc(O)cc1`) 자체도 한 행 있었고, SR-ARE와 NR-AhR에서 활성이었습니다. NR-AhR은 제가 세운 가설에 없던 결과라 이유는 확인하지 못했습니다.

코드는 `domain_analysis.py`, 전체 결과는 `domain_enrichment.csv`에 있습니다.

### 이 결과를 읽을 때 주의할 점

- 유기주석은 16개뿐이라 비율이 크게 흔들릴 수 있습니다.
- 구조 11개 그룹과 타깃 12개를 모두 비교했기 때문에, 우연히 높게 나온 조합도 섞여 있습니다. 위 표는 미리 정한 가설만 옮긴 것입니다.
- 여기서 활성은 세포 시험에서 신호가 나왔다는 뜻이고, 사람에게 독성이 있다는 증명은 아닙니다. 세포가 죽어서 생긴 신호일 수도 있습니다.
- 구조가 원인이라는 증명이 아니라 같이 나타난다는 관찰입니다.

### 구조 표시를 모델에 넣어 봤을 때

위 구조 11개를 0/1 열로 만들어 기존 특성 2,265개에 더하고(총 2,276개), 같은 방식의 scaffold 5-fold로 다시 평가했습니다. 코드는 `domain_features_cv.py`, 결과는 `domain_features_cv.csv`입니다.

12개 타깃 평균 ROC-AUC는 0.7945에서 0.7937로 거의 같았습니다. 타깃별로는 6개가 조금 올랐고 6개가 조금 내려갔으며, 12개 모두 변화량이 fold 간 표준편차보다 작아서 구조 표시를 넣어 좋아졌다고 말할 수 없습니다. SR-MMP는 5개 fold 모두 올랐지만 +0.005로 작았습니다.

Morgan 지문과 RDKit 기술자에 이미 비슷한 구조 정보가 들어 있어서 새로 더해진 정보가 없었을 수 있다고 생각합니다. 이 부분은 확인하지 못했습니다.
