# Reproducibility of a literature-derived miRNA panel for tumour-vs-normal discrimination in two independent breast-cancer cohorts

A cross-cohort benchmark of a pre-specified 15-microRNA panel on two independent Agilent miRNA microarray datasets, with explicit null comparators, calibration analysis and threshold sweep. Everything is public-data, fully reproducible from the notebook alone.

> **Scope.** This project asks whether a panel of 15 miRNAs curated from a recent breast-cancer treatment-resistance review reproduces a tumour-vs-adjacent-normal signal across two independent cohorts, *and whether it does so beyond what random 15-miRNA panels or well-known single miRNAs (miR-21, miR-210) would achieve on the same task*. It does **not** validate treatment-resistance biology, predict chemotherapy response, or claim clinical diagnostic performance. Those questions require datasets with treatment metadata and outcome, which neither GSE45666 nor GSE38167 provides.

## The question

Given how easy tumour vs adjacent-normal separation is on bulk miRNA arrays, the interesting question is not "can we hit AUC 0.98?" but:

1. Does a literature-curated panel outperform random 15-miRNA panels drawn from the same expressed pool?
2. Does the 15-marker panel add anything over the best single marker (miR-21 alone)?
3. Is the observed external AUC above the label-permutation null?
4. Are the model's probabilities calibrated across cohorts, or is a high AUC hiding a poor threshold behaviour?

The notebook answers all four.

## Datasets

| Cohort | Accession | Platform | Tumour | Adjacent-normal | Other | Role |
|---|---|---|---|---|---|---|
| Development | GSE45666 | Agilent miRNA array | 101 breast tumours | 15 | - | Model fitting, grouped 5-fold CV |
| External (locked) | GSE38167 | Independent Agilent design | 31 primary TNBC (IDC) | 23 NAT | 13 LN metastases (excluded from primary endpoint) | Locked external validation |

Both series are downloaded programmatically via `GEOparse` at run time. Sample labels are derived only from GEO titles - no manual curation.

## Method (one paragraph)

Fifteen candidate miRNAs are pre-specified from Gheorghe et al. 2025. Both cohorts are quantile-normalised in log2 space, collapsed to unique miRNA IDs by median, and converted to within-sample percentile ranks to reduce platform-level scale differences. Candidates are mapped to each platform via a family-level match (`miR-200a`, `miR-200a*`, `miR-200a-1` collapse to one median rank per family). An L2-regularised, class-balanced logistic regression on standardised rank features is trained on GSE45666 with a `StratifiedGroupKFold` grouping variable that keeps samples sharing a source identifier in the same fold, and applied unchanged to GSE38167. A patient-cluster bootstrap (2,000 draws over patient IDs) gives the external CI.

The reframed version of the notebook then adds three null benchmarks and a calibration/threshold analysis, described below.

## Results at a glance

**Cross-cohort reproducibility of the literature panel.** External ROC-AUC 0.987 (95% patient-cluster bootstrap CI 0.961-1.000, upper bound at the estimation ceiling given n = 54). Development grouped-OOF AUC 0.997 with balanced accuracy 0.985 (fragile - only 15 negatives in dev). Direction of tumour-vs-normal effect concordant across cohorts for 13 of 15 candidates; miR-34a and miR-503 flip in TNBC.

**Panel vs random-panel null (n = 1,000 random 15-miRNA draws from the shared expressed pool, candidates excluded).** Random-panel AUC median 0.881, 95th percentile 0.987. Literature panel sits at percentile **94.7**, one-sided empirical **p = 0.053**. The specific literature curation gives a modest, borderline-significant lift over an arbitrary 15-miRNA panel. The dominant driver of AUC 0.987 is the ease of tumour-vs-adjacent-normal separation itself, not the panel identity.

**Panel vs best single marker.** miR-21 alone: external AUC **0.900**. miR-210 alone: 0.872. Full 15-panel gain over miR-21 alone: **+0.087**. The panel is largely a miR-21 / miR-210 / miR-200-family story with modest additive contribution from the remaining candidates.

**Label-permutation null (n = 5,000 shuffles of external labels).** Null AUC mean 0.500, 95th percentile 0.630. Empirical p < 2 &times; 10<sup>-4</sup>. The observed AUC is not compatible with the sharp null of no signal.

**Calibration.** Brier score 0.193. Calibration slope 0.68 (ideal = 1.0), intercept 2.23 (ideal = 0.0). The model systematically under-predicts tumour probability on the external cohort due to class-prior shift between development (87% tumour) and external (57% tumour). At the default 0.5 threshold, sensitivity is 0.58 and specificity 1.00; the balanced-accuracy-optimal threshold is **0.05**, at which sensitivity rises to ~0.95 while specificity stays at 1.00. Recalibration on the target cohort - or use of percentile ranks rather than raw probabilities - is required before any threshold-based downstream use.

**One-sentence take.** A pre-specified 15-miRNA panel curated from the recent breast-cancer treatment-resistance literature reproduces the tumour-vs-adjacent-normal signal across two independent Agilent cohorts with external AUC 0.987, modestly above what random 15-miRNA panels achieve (percentile 94.7, p = 0.053), and largely driven by miR-21, miR-210, and the miR-200 family; probabilities require cohort-specific recalibration before any threshold-based use.

## Repository structure

```
mirna_crosscohort/
├── README.md                              # this file
├── mirna_crosscohort_reframed.ipynb       # main analysis notebook
├── requirements.txt                       # pinned dependencies
├── results/                               # generated on run
│   ├── GSE45666_metadata.csv
│   ├── GSE38167_metadata_all.csv
│   ├── candidate_platform_mapping.csv
│   ├── GSE45666_literature_features.csv
│   ├── GSE38167_literature_features.csv
│   ├── candidate_replication_results.csv
│   ├── development_oof_predictions.csv
│   ├── external_GSE38167_predictions.csv
│   ├── final_model_coefficients.csv
│   ├── null_random_panels.csv
│   ├── single_marker_baselines.csv
│   ├── external_threshold_sweep.csv
│   ├── results_summary.json
│   └── results_summary_extended.json
└── figures/                               # generated on run
    ├── cross_cohort_effect_replication.png
    ├── roc_external_validation.png
    ├── model_coefficients.png
    ├── null_random_panels.png
    ├── null_label_permutation.png
    ├── external_calibration.png
    └── external_threshold_sweep.png
```

## How to reproduce

The notebook was written for Kaggle Notebooks with Internet ON (so GEOparse can pull the two GEO series). It runs unchanged locally as well.

```bash
git clone https://github.com/<your-username>/mirna-crosscohort-benchmark.git
cd mirna-crosscohort-benchmark
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute mirna_crosscohort_reframed.ipynb
```

Total run time is ~5-10 minutes on a modest machine, dominated by the 1,000 random-panel refits in Section 13 and the GEO downloads in Section 2. If you want to skip the random-panel benchmark for a quick smoke test, set `N_RANDOM_PANELS = 50` in that cell.

Random seeds are fixed (`SEED = 2026`). Every intermediate table used in a figure is written to `results/` so downstream reproducibility does not depend on re-running the fits.

## Requirements

```
numpy>=1.26
pandas>=2.1
scipy>=1.11
scikit-learn>=1.4
matplotlib>=3.8
GEOparse>=2.0
jupyter>=1.0
```

Python 3.10 or later.

## What the null benchmarks are for

The single most common failure of miRNA-panel papers is claiming that a specific curated list is doing biological work when the underlying task is easy enough that most 15-marker draws would score comparably. Three null benchmarks are included so the reader can judge the panel on its own merits:

1. **Random 15-miRNA panels (Section 13).** One thousand random draws from the shared-expressed pool (excluding our candidates), each refit end-to-end on GSE45666 and evaluated on GSE38167. Produces the null AUC distribution and the empirical p-value.
2. **Single-marker baselines (Section 14).** Fifteen one-feature logistic regressions, one per candidate. Reports the best individual miRNA and the panel's marginal gain over it.
3. **Label-permutation test (Section 15).** Five thousand shuffles of the external labels with the model's predictions held fixed. Confirms the observed AUC is not compatible with the sharp null of no predictive signal.

The reader should form an opinion about the panel by looking at how these three land, not by looking at "AUC 0.987" alone.

## Calibration and threshold behaviour

Section 16 reports the Brier score, reliability curve, calibration slope and intercept on the external cohort. Section 17 sweeps the decision threshold and reports sensitivity, specificity, PPV, NPV, balanced accuracy and Youden's J across the range 0.05-0.95. This exposes the fact that the default 0.5 cutoff is not the correct operating point on the external cohort and shows what a recalibration would need to change.

## Limitations (upfront, not hidden at the end)

1. Secondary analysis of retrospective public microarray datasets from early-2010s Agilent designs. No RNA-seq or qPCR confirmation.
2. Endpoint is tumour vs adjacent-normal tissue - a near-trivial task on bulk miRNA data. Not a diagnostic or prognostic endpoint.
3. Adjacent-normal tissue is not equivalent to healthy-donor breast tissue and carries field-cancerisation effects.
4. GSE38167 is TNBC-only while GSE45666 is broader breast cancer. Cross-subtype transportability is not tested.
5. Only 15 negatives in development; the near-perfect OOF AUC is fragile.
6. Only 54 external samples; the bootstrap CI upper bound touches 1.0.
7. Two Agilent designs and historical miRBase nomenclature required probe-to-family harmonisation, which may inflate signal by combining probes.
8. Poor external calibration at threshold 0.5 - AUC is preserved but probabilities are not.
9. No treatment metadata, no clinical outcomes.
10. This is an exploratory reproducibility benchmark, not a biomarker validation study in the regulatory sense.

## What this project is a good example of

- Locking a feature panel from literature before touching performance.
- Grouped cross-validation and patient-cluster bootstrap for small, correlated samples.
- Cross-platform harmonisation via within-sample percentile ranks.
- Honest external evaluation with calibration and threshold sweep, not just AUC.
- Explicit null comparators for a task that is otherwise easy to overclaim.
- Full reproducibility from public data with a single command.

## What this project is not a good example of

- Anything about treatment response, resistance biology or clinical utility.
- Anything about miRNA biomarker discovery - every one of the 15 miRNAs was already well-known.
- Anything that could substitute for a prospective validation with real outcome data.

## Citation

If this repository was useful to you, please cite it as:

```
<Author>. Reproducibility of a literature-derived miRNA panel for tumour-vs-normal
discrimination in two independent breast-cancer cohorts. GitHub, 2026.
https://github.com/<your-username>/mirna-crosscohort-benchmark
```

Underlying datasets:

- Chen J. et al. GSE45666. Gene Expression Omnibus.
- Boissiere-Michot F. et al. GSE38167. Gene Expression Omnibus.
- Gheorghe A.-S. et al. 2025 review on miRNAs in breast-cancer treatment resistance (see notebook for full reference).

## License

MIT. See `LICENSE`.
