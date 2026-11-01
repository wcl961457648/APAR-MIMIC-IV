# APAR-MIMIC-IV — analysis code (Python)

Repository: https://github.com/wcl961457648/APAR-MIMIC-IV


Reproduction code for *"Alkaline phosphatase-to-albumin ratio and 30-day and
90-day mortality in critically ill patients with cirrhosis and sepsis"*
(cohort: MIMIC-IV v3.1, **N = 1,865**, revised 24-hour measurement window).

Everything reported in the manuscript — every table, every figure, every
number in the text — is produced by the scripts in this folder and nothing
else. There is no multiple imputation step: missingness of the model
covariates is below 2%, so the adjusted models (Models 2–4) are fitted
on complete cases (n = 1,831 of 1,865; see `09_missingness.py` and the
Missing data subsection of Material and Methods).
Model 1 is unadjusted and therefore uses all 1,865 patients (599 deaths at 30 days, 755 at 90 days); the complete-case sample has 590 and 744.

## 0. Numbering

The submitted manuscript has **6 tables and 7 figures**, numbered without an
"S" prefix. This README keeps the original script-level labels
(`TableS5_sensitivity.csv`, `12_figureS1_dag.py`, ...) because they are file
names, and gives the manuscript number of every output.

| Output | Submitted manuscript |
|---|---|
| baseline characteristics | Table 4 |
| Cox models, continuous / by quartile | Table 5 (A) / (B) |
| sensitivity analyses | Table 2 |
| timing of the index draws | Table 3 |
| discrimination, incremental value | Table 6 (a) / (b) / (c) |
| missingness | Methods, Missing data (text; no table) |
| PH tests, VIF, RMST, time-split / time-varying, interactions | Results (text; no table) |
| directed acyclic graph | Figure 5 |
| calibration | Figure 6 |
| decision-curve analysis | Figure 7 |
| flow diagram, RCS, Kaplan-Meier, subgroups | Figures 1-4 (unchanged) |

Two older labels survive inside the script docstrings and inside the SQL
comments. "Supplementary File 1" always means the cohort-extraction SQL
(`Supplementary_File1_CohortExtraction.sql` at the repository root, byte-identical
to `sql/01_cohort_extraction.sql`), and "Supplementary File 3" means the Python
analysis code, that is, this repository itself. The labels are left in place
because they sit in working code; the numbers that appear in the manuscript are
the ones in the table above.

---

## 1. What is in this folder

```
Supplementary_File3_Python/
├── README.md                  this file
├── apar_common.py             shared cohort loading, variable definitions,
│                              model specifications (M1-M4), horizons
├── 00_table1.py               Table 4   baseline characteristics
├── 01_main_models.py          Table 5    Cox models, APAR continuous and by
│                              quartile, models M1-M4, 30 d and 90 d
├── 02_discrimination.py       Table 6  C-index of single markers and of the
│                              nested models (SOFA, MELD-Na, M3, M4, each with and without APAR; shift-corrected bootstrap CIs),
│                              incremental ΔC / continuous NRI / IDI,
│                              optimism-corrected C
├── 03_calibration_dca.py      Figure 6 (calibration), Figure 7 (decision
│                              curve analysis), calibration slope and CITL
├── 04_assumptions.py          proportional-hazards tests and variance
│                              inflation factors (reported in the Results text)
├── 05_rmst_timevarying.py     restricted mean survival time, time-split
│                              intervals and time-varying-coefficient models
│                              (reported in the Results text)
├── 06_subgroups.py            Figure 4  subgroup analyses and
│                              tests for interaction
├── 07_sensitivity.py          Table 2  sensitivity analyses (complete case,
│                              3-SD trimming, competing definitions, E-value,
│                              inclusion of deaths within 24 h)
├── 08_figures.py              Figure 2  (restricted cubic splines)
│                              Figure 3  (Kaplan-Meier by APAR quartile)
├── 09_missingness.py          missingness of the candidate variables
│                              (Methods, Missing data; no table)
├── 10_timing.py               Table 3  timing of the index ALP / albumin
│                              draws (reads the cohort extract; no DB needed)
├── 11_figure1_flowchart.py    Figure 1  flow of participants
├── 12_figureS1_dag.py         Figure 5 directed acyclic graph (cohort-
│                              independent: the confounder clusters mirror the
│                              Methods and the a priori assumptions)
├── 13_calibration_optimism.py  Optimism correction for the calibration slope
│                              and calibration-in-the-large (1,000 resamples),
│                              the companion of the C-index correction in
│                              02_discrimination.py
├── 14_pairwise_logrank.py     Figure 3 pairwise log-rank tests between the
│                              APAR quartiles (all six pairs per horizon):
│                              Q1 and Q2 are indistinguishable, and the curves
│                              separate only from Q3 onwards
├── data/                      inputs you generate yourself (empty here —
│                              MIMIC-IV cannot be redistributed); see data/README.md
├── results/                   every CSV / PDF / PNG / TIFF written by the above
└── sql/
    ├── 01_cohort_extraction.sql        the script run against MIMIC-IV to
    │                                   build the analytic cohort (24-h window)
    ├── 03_sensitivity_cohorts.sql      auxiliary extracts used by Table 2
    └── 04_table_s2_timing.sql          the timestamp queries behind Table 3
```

The cohort-definition SQL is distributed twice on purpose: fully commented as
`sql/01_cohort_extraction.sql` and, at the repository root, as
`Supplementary_File1_CohortExtraction.sql`. The two are byte-for-byte
identical (same MD5).

## 2. Input data

The scripts read two spreadsheets:

| File | Content | Produced by |
|---|---|---|
| `cohort_v2.xlsx` | the analytic cohort (N = 1,865): demographics, severity scores, first measurements **within 24 h of ICU admission**, outcomes | `sql/01_cohort_extraction.sql` (byte-identical to the root-level copy) |
| `derived_v2.xlsx` | cirrhosis aetiology, AKI, ACLF grade, organ-failure components | derived from the same extract |

`09_missingness.py` additionally reads the raw extract `cohort_v2_raw.csv`
(same 1,865 patients, pre-recode column names) and `10_timing.py` reads
`cohort_v2_extras.csv` (ALP / albumin sampling times).

Both are derived from MIMIC-IV and are therefore **not redistributed here**;
access requires PhysioNet credentialing
(https://physionet.org/content/mimiciv/3.1/). See `data/README.md` for the full
list of required inputs and the query that builds each one.

### Where the scripts look for inputs

No absolute path is hard-coded. Every path is resolved by `resolve_path()` /
`resolve_outdir()` in `apar_common.py`, in this order:

1. **environment variable** — `APAR_COHORT`, `APAR_DERIVED`, `APAR_RAW`,
   `APAR_EXTRAS`, `APAR_CC_STAY`, `APAR_ALBUMIN`, `APAR_DSOFA`, `APAR_ANCHOR`,
   `APAR_BILIARY`, `APAR_COHORT_INCL24H`, `APAR_OUTDIR`;
2. **`./data/<filename>`** inside this package (simply drop the files there);
3. **`$APAR_LEGACY_DIR/<filename>`** — optional fallback pointing at a directory
   that holds the same file names (e.g. your own copy of the extraction folder).

Command-line flags always win. So either of these works:

```bash
# (a) environment variables
export APAR_COHORT=/path/to/cohort_v2.xlsx
export APAR_DERIVED=/path/to/derived_v2.xlsx
export APAR_OUTDIR=./results
python 01_main_models.py

# (b) explicit flags
python 01_main_models.py --cohort /path/to/cohort_v2.xlsx \
                         --derived /path/to/derived_v2.xlsx \
                         --outdir ./results
```

Every script accepts `--cohort`, `--derived`, `--outdir` and `--n-boot`
(plus the file-specific flags listed in `data/README.md`).

A file that cannot be found is skipped with a `[skip]` message instead of
aborting — only the corresponding sensitivity analysis is dropped from Table 2.

## 3. Software

* **Cohort extraction** — PostgreSQL 11.19 against MIMIC-IV v3.1.
* **All analyses** — Python 3.13.14. Versions used to generate every file in
  `results/`: NumPy 2.5.3, pandas 2.3.3, SciPy 1.18.1, lifelines 0.30.3,
  statsmodels 0.15.0, Matplotlib 3.11.1, Pillow 12.3.0 (flow chart and DAG
  figures), openpyxl 3.1.5 (Excel input).
* **Roles** — lifelines: Cox models, Kaplan–Meier, Schoenfeld residuals,
  Harrell's C-index; statsmodels: collinearity diagnostics (VIF); Matplotlib:
  all figures. **scikit-survival is not used** — Harrell's C-index comes from
  `lifelines.utils.concordance_index` and `CoxPHFitter.concordance_index_`.
* **Confidence intervals for hazard ratios** — Wald intervals on the log
  hazard ratio scale, `exp(beta ± 1.96 x SE(beta))` (see `01_main_models.py`
  lines 64-66 and 87-89). lifelines offers no profile-likelihood option and
  none is used. By contrast, the bootstrap CIs (ΔC-index, NRI, IDI,
  calibration slope, calibration-in-the-large) are percentile intervals.
* **Decision curve analysis** — implemented directly from the net-benefit
  definition in `03_calibration_dca.py::net_benefit` (plain NumPy; no external
  DCA package). The outcome is dichotomised at each horizon
  (`(T <= t_eval) & (E == 1)`); inverse-probability-of-censoring weighting is
  unnecessary because no patient was lost to follow-up before 30 or 90 days.
* No imputation software is used (complete-case analysis; see above).

```bash
pip install numpy pandas scipy lifelines statsmodels matplotlib openpyxl pillow
```

## 4. Running the code in order

```bash
python 00_table1.py            # seconds
python 01_main_models.py       # seconds
python 02_discrimination.py    # slow: 1,000 bootstrap resamples (~20 min)
python 03_calibration_dca.py   # ~10 min
python 04_assumptions.py       # seconds
python 05_rmst_timevarying.py  # ~1 min
python 06_subgroups.py         # ~1 min
python 07_sensitivity.py       # ~2 min
python 08_figures.py           # ~1 min
python 09_missingness.py       # seconds
python 10_timing.py            # seconds (reads cohort_v2_extras.csv)
python 11_figure1_flowchart.py # seconds
python 12_figureS1_dag.py      # seconds
python 13_calibration_optimism.py # ~6 min
python 14_pairwise_logrank.py    # seconds
```

For a fast smoke test, add `--n-boot 50`; note that the optimism correction in
`02_discrimination.py` needs the full 1,000 resamples to be stable — with 50
resamples the mean optimism is smaller than its own Monte-Carlo error and the
sign of the correction is unreliable.

## 5. Conventions that matter for reproduction

These points silently change the numbers if you deviate from them:

1. **24-hour window.** ALP, albumin and every laboratory covariate are the
   *first* measurement with `icu_intime <= charttime < icu_intime + 24 h`.
   Severity scores (SOFA, APS III, OASIS, SAPS II, GCS) are the official
   first-day derived values, whose window runs from 6 h before to 24 h
   after ICU admission rather than from 0 to 24 h.
   The Glasgow Coma Scale used for the brain-failure component of ACLF is the
   first-day minimum (`first_day_gcs`), so that component shares this wider
   window; restricting it to 0–24 h would give 155 brain failures
   instead of 157. The cohort also carries `gcs_whole`, the minimum over the
   entire ICU stay; it is used **only** by row 15 of Table 2 and by nothing
   else in the primary analysis.
2. **Units and APAR scale.** APAR = ALP (IU/L) ÷ albumin (**g/L**). MIMIC-IV
   stores albumin as g/dL, so the stored value is multiplied by 10 before the
   division (see `apar_from()` in `apar_common.py`); albumin itself is likewise
   reported in g/L everywhere in the manuscript and tables. Because the factor
   is a pure scale change, every P value, quartile hazard ratio, per-1-SD
   hazard ratio, C-index, NRI, IDI and RMST is unaffected — only the HR per one
   APAR unit and the numbers printed on the APAR scale move. The two analytes
   need not be drawn at the same time; each contributes its own first value in
   the window (88.2% of patients had both drawn at the same time; see Table 3).
3. **No outlier deletion in the primary analysis.** Patients with APAR above
   mean + 3 SD (threshold 18.89, n = 36) are retained; the exclusion is
   reported as sensitivity analysis 3 only.
4. **INR** is the laboratory-reported international normalised ratio
   (MIMIC-IV item 51237).  It is *not* approximated as prothrombin time (s)
   / 11.0: that approximation holds only if the international sensitivity
   index is 1 and 11.0 s is exactly the laboratory's normal mean prothrombin
   time, neither of which can be verified in MIMIC-IV.  The two agree exactly
   in only 8.4% of patients (Pearson r = 0.9722) and the normal mean PT
   implied by the reported INR has a median of 10.889 s; PT and INR come from
   the same charttime in 1,831 of 1,833 patients (99.89%), so the difference
   is not a timing artefact.  Using the reported INR also makes MELD-Na
   consistent with the ACLF coagulation criterion, which has always used item
   51237.  The PT / 11.0 s variant is retained as sensitivity analysis 16 in
   Table 2 (`07_sensitivity.py`).
5. **MELD / MELD-Na — full 2016 OPTN bounding rules.** All four OPTN bounds are
   applied (OPTN Policy 9, effective January 2016):
   ```
   Cr_c    = min( max(creatinine, 1.0), 4.0 )     # floor 1.0, then ceiling 4.0
   Bili_c  = max(bilirubin, 1.0)                  # floor 1.0
   INR_c   = max(reported INR, 1.0)               # floor 1.0
   MELD(i) = max(0.957 ln(Cr_c) + 0.378 ln(Bili_c) + 1.120 ln(INR_c) + 0.643, 0) × 10
   Na_c    = min(137, max(125, sodium))           # inclusive at both ends
                                                  # ≤ 125 → 125, ≥ 137 → 137
   MELD-Na = min( max( MELD(i) + [MELD(i) > 11]
                       × (1.32(137 − Na_c) − 0.033·MELD(i)(137 − Na_c)), 0 ), 40 )
   ```
   The sodium clamp is continuous, so no observed value can fall into an
   undefined interval. `apar_common._self_test_meld_bounds()` asserts every bound
   on import (Na 120→125, 137→137, 140→137; Cr 0.2→1.0, 18.1→4.0; Bili 0.1→1.0;
   INR 0.8→1.0; and the MELD > 11 gate), and `meld_na_from()` in `apar_common.py`
   is the single implementation (`07_sensitivity.py` calls it as well). The
   identical rules are written in `sql/01_cohort_extraction.sql` (STEP 8,
   `t_meld`).

   **One documented deviation.** OPTN rounds MELD(i) to an integer for
   allocation; here it is kept continuous because it is used as a continuous
   covariate and as a subgrouping variable. Rounding changes no bound and is
   not applied by most prognostic-research implementations of the score.
6. **C-index confidence intervals.** Bootstrap C of a model fitted to a
   bootstrap resample is biased downwards (duplicated patients reduce the
   proportion of comparable pairs), so the bootstrap distribution is recentered
   on the point estimate before the percentiles are taken:
   `CI = percentiles(bs − (median(bs) − C_observed))`.  The same recentering is
   used for the paired ΔC-index.
7. **Bootstrap optimism.** The bootstrap frame *must* be
   `reset_index(drop=True)` before fitting: lifelines' `concordance_index_`
   collapses towards 0.5 when the fitted frame carries duplicated index labels.
   Harrell's estimator is used: optimism = mean over B resamples of
   [C(model_b on its own resample) − C(model_b on the original sample)],
   corrected C = apparent C − optimism.
8. **Administrative censoring.** Follow-up is censored at 30 and 90 days after
   ICU admission; `dod` is a date, so survival time is counted in whole days.
   "Events" for the NRI at horizon *t* are therefore deaths with
   `T ≤ t` and `E = 1`; everyone else is a non-event.

   The same convention governs the number-at-risk row printed under
   each Kaplan-Meier panel: it counts the patients with `T ≥ t`, so at
   the horizon it still contains the deaths recorded on that exact day.
   Reconstructing the count as `n - events` therefore understates it by
   1-3 at 30 days (Q1 to Q4: 1, 2, 3 and 1 deaths on day 30) and by 1 at
   90 days (Q3 only), while no patient is censored before either horizon.
   `Figure3_km_summary.csv` reports the two components as `deaths_at_tau`
   and `censored_before_tau`, so that
   `at_risk_at_tau = n - events + deaths_at_tau - censored_before_tau`
   can be checked without re-deriving the survival times.

9. **Outcome column names.** `death_within_icu_30days` and
   `death_within_hosp_90days` are legacy labels inherited from the earlier
   N = 2,209 extraction; in this cohort they are **not** setting-specific.
   Both carry all-cause mortality derived from `patients.dod` (SQL STEP 7) and
   are row-for-row identical to `death_30d` / `death_90d`; `load_cohort()`
   asserts that identity at run time. Because `dod` merges hospital and state
   death-registry records, deaths occurring after discharge are captured (617
   patients here have dod after the index discharge, 239 within 90 days).
   In-hospital (505) and ICU (326) mortality are reported in Table 4 as
   supportive outcomes and are not modelled anywhere.

## 6. Random seeds

All resampling uses fixed seeds (`20260910` for optimism, `11` for NRI/IDI,
`42` for calibration, `5` for ΔC-index between models and `3` for the C-index of
single predictors), so re-running the scripts reproduces the reported values to
the printed precision. Bootstrap-based intervals (ΔC-index, NRI, IDI,
calibration slope) differ by <0.002 between runs with different seeds.

The number of resamples is not a single number but is fixed per purpose, and is
reported as such in the Methods:

| Purpose | Resamples | Where |
|---|---|---|
| C-index of single predictors; ΔC-index, NRI, IDI (paired, recentered percentile CIs) | 300 | `02_discrimination.py` (`--n-boot-nri`, `--n-boot-cindex`) |
| Calibration slope and calibration-in-the-large (plain percentile CIs) | 500 | `03_calibration_dca.py` (`--n-boot-cal`) |
| Optimism correction of the C-index (a mean over resamples, which enters the point estimate) | 1,000 | `02_discrimination.py` (`--n-boot-opt`) |

The decision curve analysis uses no resampling at all: net benefit is computed
directly on the original sample (`03_calibration_dca.py::net_benefit`).

## 7. Output ↔ manuscript crosswalk

| Output file | Manuscript item |
|---|---|
| `results/Figure1_flowchart.*` | Figure 1 (flow diagram); the exported image carries no in-figure title or footnote — both live in the manuscript legend |
| `results/Table1_baseline.csv` | Table 4 |
| `results/Table2_continuous_APAR.csv` | Table 5 (A) |
| `results/Table3_quartiles_APAR.csv` | Table 5 (B) |
| `results/Table3_quartile_counts_complete_cases.csv` | Table 5 (B): per-quartile patients and deaths in the 1,831 complete cases in which Models 2–4 are fitted |
| `results/Table4_Cindex.csv` | Table 6 (a): C-index (95% CI) of every single marker and of the nested models |
| `results/Table4_incremental_NRI_IDI.csv` | Table 6 (b) |
| `results/Table4_optimism_corrected_Cindex.csv` | Table 6 (c): optimism-corrected C-index, 16 rows (SOFA, MELD-Na, Model 3 and Model 4, each with and without APAR, at 30 and 90 days) |
| `results/e_value.csv`, `results/mortality_by_quartile.csv` | text, Results |
| `results/Figure2_RCS.*`, `Figure2_rcs_nonlinearity.csv` | Figure 2 |
| `results/Figure3_KM.*`, `Figure3_km_summary.csv` | Figure 3; the summary file holds one row per horizon and quartile (n, events, Kaplan-Meier survival at the horizon, number at risk at the horizon, and the two components `deaths_at_tau` and `censored_before_tau` that reconcile the risk count with the death count) |
| `results/Figure4_subgroup_*.{pdf,png,tiff}`, `Figure4_subgroup_estimates.csv` | Figure 4 |
| `results/TableS1_missingness.csv` | Methods, Missing data (text; no table) |
| `results/TableS2_timing.csv` | Table 3 (sampling-time distribution) |
| `results/TableS3_proportional_hazards.csv` | Results text (proportional-hazards) |
| `results/TableS4_vif.csv` | Results text (variance inflation factors) |
| `results/TableS5_sensitivity.csv` | Table 2 |
| `results/TableS6_interactions.csv` | Results text (subgroup interactions; Figure 4 is the figure) |
| `results/TableS7_rmst_*.csv` | Results text (RMST differences) |
| `results/TableS8_timesplit_*.csv`, `TableS8_tvc_hrs.csv`, `TableS8_tvc_tests.csv` | Results text (time-split and time-varying models; `TableS8_tvc_tests.csv`
holds the P for interaction with log(time) quoted in the Results; it is a Wald test on the `apar_z:g` term) |
| `results/FigureS1_DAG.*` | Figure 5 (directed acyclic graph); the exported image carries no in-figure title or footnote — the Model 3 adjustment set and the >15%-missingness rule are repeated in the manuscript legend |
| `results/calibration_slope_optimism.csv` | Figure 6 and text, Results: optimism-corrected calibration slope and calibration-in-the-large |
| `results/FigureS2_calibration.*`, `calibration_slope_CITL.csv` | Figure 6 (the panels keep their 30/90-day axis titles; only the figure-level title was dropped) |
| `results/FigureS3_DCA.*` | Figure 7 (idem) |

Figure 4 is written as one wide figure: the 30-day panel (A) and the
90-day panel (B) sit side by side and share a single subgroup column,
because the subgroup membership and the group sizes are identical at the
two horizons and only the event counts and the HRs differ.  The figure is
drawn natively 180 mm wide (2126 px at 300 dpi) and is never resampled, so
the 7-9.5 pt type is the type that prints; it is about 180 x 154 mm.  Two
labels wrap onto a second line ("Alcohol-related aetiology", "Other
(non-alcoholic) aetiology") rather than being shortened, so the figure text
stays identical to Figure4_subgroup_estimates.csv.  The estimates carry an
"HR (95% CI)" column header and the P values a "P" header, which is why the
values are printed without a "p =" prefix; three decimals are used and the
value is replaced by the bound "<0.001" below that threshold, because two
decimals would render every P between 0.001 and 0.005 as "0.00".  The x
axis is ticked at 1.0 and 1.5 only: three ticks in log coordinates left the
1.0 and 1.2 labels 12 pt apart on a 17-mm axis and merged into "1.01.2".

## 8. Table 2 — the rows

Row numbers below are the labels written by `07_sensitivity.py`, which are also
the numbers used in Table 2 of the manuscript: 1 is the first sensitivity row.
The two rows above it — the primary Model 3 analysis and the complete-case
analysis — are unnumbered in both places, and letters mark sub-rows (4a / 4b,
7b / 7c, 13a–13c).

| Row | Sensitivity analysis | n | Q4 vs Q1, 30 d | Q4 vs Q1, 90 d |
|---|---|---|---|---|
| — | Primary analysis (Model 3) | 1,831 | 1.685 (1.344–2.112) | 1.851 (1.512–2.266) |
| — | Complete-case analysis | 1,831 | 1.724 (1.375–2.163) | 1.884 (1.539–2.308) |
| 1 | APAR > mean + 3 SD excluded (threshold 18.89, n = 36) | 1,795 | 1.680 (1.335–2.113) | 1.790 (1.458–2.198) |
| 2 | APAR winsorised at the 1st/99th percentile | 1,831 | 1.685 (1.344–2.112) | 1.851 (1.512–2.266) |
| 3 | log2(APAR), HR per doubling | 1,831 | 1.248 (1.153–1.351) | 1.281 (1.192–1.377) |
| 4a | ALP alone, per 1 SD (SD = 103.5 IU/L) | 1,831 | 1.171 (1.104–1.241) | 1.195 (1.134–1.259) |
| 4b | Albumin alone, per 1 SD (SD = 6.53 g/L) | 1,831 | 0.925 (0.857–0.999) | 0.913 (0.851–0.979) |
| 5 | Excluding albumin infusion before the index draw (164 of the 1,831 complete cases excluded; 167 cohort patients carried the flag) | 1,667 | 1.843 (1.435–2.367) | 2.048 (1.642–2.554) |
| 6 | Excluding biliary obstruction / cholangitis (150 of the 1,831 complete cases excluded; 152 cohort patients carried the flag) | 1,681 | 1.732 (1.367–2.195) | 1.853 (1.500–2.289) |
| 7 | Sepsis redefined by a baseline-corrected SOFA rise ≥ 2 | 1,459 | 1.874 (1.431–2.453) | 2.143 (1.683–2.730) |
| 7b | Idem, only the 379 patients meeting the criterion | 379 | 1.986 (1.164–3.389) | 1.951 (1.218–3.124) |
| 7c | P for interaction between baseline availability and APAR (not an HR) | 1,831 | 0.510 | 0.232 |
| 8 | E-value (point estimate / CI limit) | — | 2.22 / 1.76 | 2.43 / 2.00 |
| 9 | Excluding admissions from 2020 onwards | 1,532 | 1.903 (1.480–2.447) | 1.942 (1.554–2.427) |
| 10 | Further adjusted for the white-cell count (Model 4; lactate is not entered into any model because it is missing in 18.5% of the cohort) | 1,831 | 1.592 (1.266–2.002) | 1.765 (1.439–2.166) |
| 11 | Deaths within 24 h of ICU admission retained (N = 1,927), **per 1 unit** | 1,891 | 1.036 (1.025–1.048) | 1.041 (1.030–1.051) |
| 12 | Time origin moved from ICU admission to the index draw (immortal-time analysis) | 1,831 | 1.677 (1.338–2.103) | 1.846 (1.508–2.260) |
| 13a | Index draw **before** sepsis onset, per 1 SD | 933 | 1.173 (1.073–1.283) | 1.221 (1.127–1.322) |
| 13b | Index draw **after** sepsis onset, per 1 SD | 898 | 1.188 (1.100–1.283) | 1.202 (1.122–1.288) |
| 13c | P for interaction between draw timing and APAR (not an HR) | 1,831 | 0.660 | 0.988 |
| 14 | INR approximated as prothrombin time / 11.0 s instead of the laboratory-reported value (item 51237) | 1,831 | 1.679 (1.339–2.104) | 1.845 (1.507–2.259) |
| 15 | Brain-failure criterion of ACLF taken over the **whole ICU stay** instead of the first 24 h | 1,831 | 1.687 (1.345–2.114) | 1.851 (1.512–2.266) |

Row 8 (E-value) is computed on the risk-ratio scale, not by substituting the
hazard ratio into the formula directly.  Mortality by the end of follow-up was
32.1% at 30 days and 40.5% at 90 days, far above the 15% threshold below which
a hazard ratio may be read as an approximate risk ratio, so the hazard ratio is
first converted to an approximate risk ratio (VanderWeele, Biometrics 2020) and
only then is E = RR + sqrt(RR * (RR - 1)) applied.  This reproduces R's
`EValue::evalues.HR(rare = FALSE)`.  Every script calls the single shared
implementation in `apar_common.e_value()`.  An earlier independent copy of the
formula in `01_main_models.py` used the rare-outcome shortcut and returned
2.76 / 3.11 instead of 2.22 / 2.43; those figures must not be used.

Row 5 definition: `mimiciv_icu.inputevents`, albumin itemids 220574 /
220861–220864 / 227456, infusion started **before** the index draw.

Row 11 is reported **per 1 unit** of APAR rather than as Q4 vs Q1, because every one
the additionally retained patients dies inside the first 24 h of follow-up; a
quartile contrast would be dominated by those immediate events.  The per-unit scale is
also the right one here: the retained cohort has a wider APAR distribution (SD 4.98
versus 4.71), so a per-1-SD estimate would not be comparable with the primary one.
The per-1-SD reading of the *same* fitted model is carried in the `Note` column of
both row-11 records, so the values quoted in the Table 2 footnote (1.195 at 30 days
and 1.220 at 90 days) can be reproduced from the supplied code.

Rows 13a–13c answer the review question about the overlap between the index
draw and sepsis onset: the draw precedes onset in 943 of 1,865 patients
(50.6%, see `10_timing.py` Panel E).  Both strata are scaled by the SD of the
whole cohort (4.706), so 13a, 13b and the primary per-1-SD estimate share one
scale; the interaction P of 0.658 (30 d) and 0.983 (90 d) shows no effect
modification by the timing of the draw.

Row 14 repeats the primary comparison with the INR approximated as
prothrombin time / 11.0 s; the two INR definitions correlate at r = 0.9722 and
give near-identical estimates, so nothing turns on how the INR is defined.

Row 15 answers the review comment about the brain-failure component of ACLF.
The primary analysis uses the minimum Glasgow Coma Scale within the first 24 h
of ICU admission, the same window as the other five components; row 15 widens
that single component to the lowest score recorded at any time in the ICU stay.
413 patients meet the criterion under the wider window and 157 under the 24-h
window, and the estimates move from 1.680 (1.340–2.107) to 1.679
(1.339–2.105) at 30 days and from 1.848 (1.509–2.262) to
1.844 (1.506–2.257) at 90 days, i.e. nowhere.  The wider window is
reported as a sensitivity analysis rather than the primary one because whether
a low score is ever recorded depends on how long the patient was observed as
well as on how ill the patient was: the 256 patients identified only under the
wider window had a median ICU stay of 10.1 days, so survival time itself
determines who is classified as having brain failure.

## 9. Reproducing the auxiliary inputs

| File | Script |
|---|---|
| `cc_stay.csv` (complete-case flag) | `sql/03_sensitivity_cohorts.sql` |
| `albumin_before_draw.csv` | `sql/03_sensitivity_cohorts.sql` |
| `anchor_year_group.csv` (calendar period) | `sql/03_sensitivity_cohorts.sql` |
| `biliary.csv` (biliary obstruction / cholangitis) | `sql/03_sensitivity_cohorts.sql` |
| `dsofa_cohort.csv` (baseline-corrected SOFA) | `sql/03_sensitivity_cohorts.sql` |
| `cohort_v2_extras.csv` (ALP / albumin sampling times) | `sql/04_table_s2_timing.sql` |
