/* ============================================================================
   REPRODUCIBLE COHORT DEFINITION AND EXTRACTION CODE  (revised 24-hour window)
   ----------------------------------------------------------------------------
   Manuscript : "Association of the alkaline phosphatase-to-albumin ratio with
                 30- and 90-day all-cause mortality in critically ill patients
                 with cirrhosis and sepsis"
   Target     : the analytic cohort of N = 1,865 patients reported in the paper
   Database   : MIMIC-IV v3.1 (PhysioNet), PostgreSQL
   Schemas    : mimiciv_hosp, mimiciv_icu, mimiciv_derived
   Runtime    : ~10-15 min on a local PostgreSQL installation
   Submitted  : Supplementary File 1 (the identical file is also shipped
                inside Supplementary File 3 as sql/01_cohort_extraction.sql).

   WHAT THIS SCRIPT DOES
   ---------------------
   It reproduces, from the raw MIMIC-IV tables, exactly the cohort used for
   every analysis in the manuscript: a single SQL file, no external input,
   ending in a table `apar_cohort_1865` with one row per patient and a final
   COUNT that must read 1,865.

   THE ITEMS ASKED FOR BY REVIEWER #1
   ----------------------------------
   (1) exact ICD codes for cirrhosis ................ STEP 1 (8 codes listed)
   (2) suspected infection window .................. STEP 3 (Appendix A)
   (3) SOFA calculation method ..................... STEP 3 (Appendix B)
   (4) baseline SOFA definition .................... STEP 3 (Appendix B)
   (5) full MIMIC-IV extraction algorithm .......... STEPS 1-9 (this file)

   WHAT CHANGED RELATIVE TO THE PREVIOUS VERSION
   ---------------------------------------------
   * THE 24-HOUR RULE.  The exposure (ALP and albumin) and every laboratory
     covariate are now the FIRST measurement taken in the 24 hours that follow
     ICU admission (charttime >= icu_intime AND charttime < icu_intime + 24 h).
     Previously the laboratory window was open-ended ("first value after ICU
     admission").  Severity scores were taken from the official first-day derived
     tables, whose window runs from 6 h before to 24 h after ICU
     admission.
   * APAR  = first ALP (IU/L) in 24 h  /  first albumin (g/L) in 24 h.
     MIMIC-IV records albumin in g/dL, so the stored value is multiplied by 10
     before the division (see `t_with_apar` below); this is the convention used
     in the published APAR literature and throughout the manuscript, where
     albumin itself is also reported in g/L.  Because it is a pure scale
     factor, every P value, quartile hazard ratio, per-SD hazard ratio,
     C-index, NRI, IDI and RMST is unaffected.
     The two analytes need NOT be drawn at the same time; each contributes its
     own first value inside the window.  The resulting sampling-time
     distribution is reported in Supplementary Table S2 and can be reproduced
     with Appendix C below.
   * ORDER OF THE EXCLUSIONS.  Patients without either an ALP or an albumin
     measurement inside the 24 h window are excluded first (3,068 -> 1,927)
     and patients dying within 24 h of ICU admission are excluded second
     (1,927 -> 1,865).
   * NO OUTLIER EXCLUSION IN THE PRIMARY ANALYSIS.  Patients with APAR above
     mean + 3 SD (n = 36, threshold 18.89) are RETAINED in the primary
     analysis; their exclusion is reported as sensitivity analysis 1 only
     (Appendix D reproduces that variant).

   HOW TO RUN
   ----------
     psql -d mimic -f cohort_extraction.sql
   (requires PhysioNet credentialed access to MIMIC-IV v3.1)

   CONVENTIONS
   -----------
   * Index ICU stay = the earliest ICU stay in which Sepsis-3 criteria were
     met; one row per patient (independence of observations).
   * "First measurement in 24 h" = the earliest charttime >= icu_intime and
     < icu_intime + 24 h, ordered by (charttime, labevent_id) so that the
     ordering is deterministic.
   * All timestamps are shifted into the future by the MIMIC-IV date-shifting
     procedure; intervals (24 h, 30 days, 90 days) are unaffected.
   ============================================================================ */

SET search_path TO mimiciv_hosp, mimiciv_icu, mimiciv_derived, public;

/* ============================================================================
   STEP 1. CIRRHOSIS - ICD-9-CM / ICD-10-CM codes (all diagnosis positions)
   ----------------------------------------------------------------------------
   Cirrhosis was ascertained from ANY diagnosis position of the index hospital
   admission (not only the principal diagnosis). The eight codes below are the
   complete set used; no other cirrhosis inclusion code was accepted.
   (K7032 is not an inclusion code: it represents alcoholic fibrosis of
   liver, not established cirrhosis, and is therefore excluded from STEP 1;
   it is used in STEP 8 solely as an aetiology marker for alcohol-related
   cirrhosis.)

     Code    Version   Meaning
     ------  --------  ------------------------------------------------------
     5712    ICD-9     Alcoholic cirrhosis of liver
     5715    ICD-9     Cirrhosis of liver, unspecified
     5716    ICD-9     Biliary cirrhosis
     K7030   ICD-10    Alcoholic cirrhosis of liver without ascites
     K7031   ICD-10    Alcoholic cirrhosis of liver with ascites
     K743    ICD-10    Primary biliary cirrhosis / cholangitis
     K7460   ICD-10    Unspecified cirrhosis of liver
     K7469   ICD-10    Other and unspecified cirrhosis of liver

   Note on the version filter.  The two code sets do not overlap, and this
   was verified against the database rather than assumed from code length:
     * 0 rows anywhere in diagnoses_icd carry an ICD-10 code on an
       icd_version = 9 row, and 0 rows carry an ICD-9 code on an
       icd_version = 10 row;
     * 0 icd_version = 9 rows begin with 'K' (in MIMIC-IV the only
       alphabetic ICD-9-CM codes are E codes and V codes).
   The version filter is nevertheless applied below for transparency.  Codes
   are stored both with and without a decimal point in MIMIC-IV, hence
   REPLACE(icd_code, '.', '').

   Empirical yield of each code (distinct hospital admissions over the whole
   of MIMIC-IV v3.1; these are database-wide counts, not the analytic cohort):

     Retained   5712 4,321 | 5715 6,100 | 5716   360
                K7030 3,021 | K7031 3,196 | K743   376
                K7460 2,876 | K7469 3,305
     Excluded   K740   166 | K741     3 | K742     3
                K744    16 | K745    31
                (K717, toxic liver disease with fibrosis, 14 admissions, is
                 not a cirrhosis code and was never part of this set)

   Asymmetry across ICD versions: ICD-9 5716 "biliary cirrhosis" does not
   distinguish primary from secondary biliary cirrhosis, whereas ICD-10 splits
   them into K743 (primary, retained) and K744/K745 (secondary, excluded).
   The ICD-9 side is therefore less specific by construction; 360 admissions
   database-wide are affected and this is not corrected for.

   Specificity is deliberately favoured over sensitivity: K740-K742 (fibrosis
   without established cirrhosis) and K744/K745 (secondary biliary cirrhosis,
   a consequence of biliary obstruction) are excluded, the latter being
   addressed instead by the prespecified sensitivity analysis that removes
   patients with biliary disease.  The imperfect sensitivity that remains is
   acknowledged as a limitation in the manuscript (Section 2.2) and is NOT
   examined by any widened-code sensitivity analysis.
   ============================================================================ */
DROP TABLE IF EXISTS t_cirrhosis;
CREATE TEMP TABLE t_cirrhosis AS
SELECT DISTINCT d.subject_id, d.hadm_id
FROM mimiciv_hosp.diagnoses_icd d
WHERE (d.icd_version = 9  AND REPLACE(UPPER(d.icd_code), '.', '')
           IN ('5712','5715','5716'))
   OR (d.icd_version = 10 AND REPLACE(UPPER(d.icd_code), '.', '')
           IN ('K7030','K7031','K743','K7460','K7469'));

-- one row per cirrhosis admission
SELECT COUNT(*) AS n_cirrhosis_admissions FROM t_cirrhosis;


/* ============================================================================
   STEP 2. STUDY POPULATION - index ICU stay, adults
   ----------------------------------------------------------------------------
   Inclusion: (i) age >= 18 years at admission;
              (ii) ONE stay per patient - the index stay is the earliest ICU stay
                   in which Sepsis-3 criteria were met.
   Age is taken from mimiciv_derived.age (anchor-age based, the standard
   MIMIC-IV age derivation).

   ORDER OF OPERATIONS (important): the "one stay per patient" restriction is
   applied to the set of stays that already satisfy the Sepsis-3 criterion
   (STEP 3), not before it.  A patient whose first ICU stay is sepsis-free but
   whose second stay meets Sepsis-3 would otherwise be dropped, or - worse -
   would contribute a stay that is not the one in which sepsis occurred.
   ============================================================================ */
DROP TABLE IF EXISTS t_icu;
CREATE TEMP TABLE t_icu AS
SELECT icu.subject_id,
       icu.hadm_id,
       icu.stay_id,
       icu.admittime,
       icu.dischtime,
       icu.icu_intime,
       icu.icu_outtime,
       icu.los_icu,
       icu.gender,
       icu.hospital_expire_flag,
       a.age
FROM mimiciv_derived.icustay_detail icu
JOIN mimiciv_derived.age a ON a.hadm_id = icu.hadm_id
WHERE a.age >= 18
  AND icu.hadm_id IN (SELECT hadm_id FROM t_cirrhosis);

SELECT COUNT(*) AS n_adult_cirrhosis_icu_stays FROM t_icu;


/* ============================================================================
   STEP 3. SEPSIS-3
   ----------------------------------------------------------------------------
   Sepsis was defined by Sepsis-3 criteria: suspected infection plus an
   increase of >= 2 SOFA points attributable to the infection.

   We use the official MIMIC-IV derived table `mimiciv_derived.sepsis3`
   (Guan et al., "Sepsis-3 definition in MIMIC-IV", the implementation
   distributed with the MIMIC Code Repository), which applies exactly those
   criteria. Its columns are used as follows:

     antibiotic_time          first systemic antibacterial agent
     culture_time             first microbiological specimen obtained
     suspected_infection_time anchor of the suspected-infection window
     sofa_time                time at which the SOFA criterion was met (t0)
     sofa_score               total SOFA at sofa_time, exported as
                              `sofa_at_onset` ("SOFA at sepsis onset").  It
                              implements the Sepsis-3 criterion (Methods 2.2)
                              and is not used as a covariate: the covariate
                              in every model is the first-day SOFA `fds.sofa`
                              added in STEP 8 (Methods 2.4, Table S10).  Same
                              six-component rule but a different window - the
                              two columns differ in 1,645 of 1,865 patients
                              (median 5 vs 9).
     respiration ... renal    the six SOFA sub-scores

   Using the official derived table (rather than re-implementing Sepsis-3)
   is deliberate: it makes the cohort definition reproducible by any
   credentialed MIMIC-IV user without re-coding the criteria.

   The two definitions that the reviewer specifically asked about are
   reproduced in full in APPENDIX A (suspected infection window) and
   APPENDIX B (SOFA calculation and baseline SOFA).
   ============================================================================ */
DROP TABLE IF EXISTS t_sepsis;
CREATE TEMP TABLE t_sepsis AS
SELECT * FROM (
    SELECT i.*,
           s3.antibiotic_time,
           s3.culture_time,
           s3.suspected_infection_time,
           s3.sofa_time          AS sepsis_onset_time,   -- t0
           s3.sofa_score         AS sofa_at_onset,
           s3.respiration, s3.coagulation, s3.liver,
           s3.cardiovascular, s3.cns, s3.renal,
           ROW_NUMBER() OVER (PARTITION BY i.subject_id
                              ORDER BY i.admittime, i.icu_intime) AS rn_first
    FROM t_icu i
    JOIN mimiciv_derived.sepsis3 s3
      ON s3.stay_id = i.stay_id
     AND s3.sepsis3 = TRUE
) t
WHERE rn_first = 1;

-- expected: 3,068 (adults with cirrhosis and Sepsis-3 during their index ICU stay)
SELECT COUNT(*) AS n_cirrhosis_sepsis3 FROM t_sepsis;


/* ============================================================================
   STEP 4. EXPOSURE AND COVARIATES - FIRST VALUE WITHIN 24 h OF ICU ADMISSION
   ----------------------------------------------------------------------------
   The window is  [icu_intime, icu_intime + 24 h)  for every laboratory item.

     itemid 50863 = Alkaline Phosphatase, IU/L   (blood)   -> exposure numerator
     itemid 50862 = Albumin, g/dL                (blood)   -> stored denominator
     (itemids 53085 / 53086 are the same analytes in urine and are not used)

   The same rule is applied to every other laboratory covariate (creatinine,
   bilirubin, sodium, PT, INR, lactate, haemoglobin, platelets, white cells,
   ...), so that exposure and confounders are measured over the identical time
   window.  Severity of illness is taken from the official first-day derived
   tables, whose window runs from 6 h before to 24 h after ICU
   admission.

   APAR = ALP (IU/L) / albumin (g/L), i.e. ALP / (albumin_g_dL x 10).
   ============================================================================ */
DROP TABLE IF EXISTS t_lab;
CREATE TEMP TABLE t_lab AS
WITH ranked AS (
    SELECT b.stay_id,
           le.itemid,
           le.valuenum,
           le.charttime,
           ROW_NUMBER() OVER (PARTITION BY b.stay_id, le.itemid
                              ORDER BY le.charttime, le.labevent_id) AS rn
    FROM t_sepsis b
    JOIN mimiciv_hosp.labevents le
      ON le.hadm_id = b.hadm_id
     AND le.itemid IN (
            50862,  -- albumin, g/dL
            50863,  -- alkaline phosphatase, IU/L
            50861, 50878, 50885, 50883,  -- ALT, AST, total / direct bilirubin
            50912, 50983,                -- creatinine, sodium
            50813,                       -- lactate
            51274, 51237,                -- PT, INR
            51222, 51265, 51301, 51279,  -- Hb, platelets, WBC, RBC
            52075, 51133,                -- absolute neutrophil / lymphocyte count
            50852, 50931,                -- HbA1c, glucose
            50820, 50818, 50821, 50882,  -- pH, pCO2, pO2, bicarbonate
            51006, 50963,                -- urea nitrogen, NT-proBNP
            50907, 50904, 50905, 51000   -- cholesterol (total/HDL/LDL), triglycerides
         )
     AND le.valuenum IS NOT NULL
     AND le.charttime >= b.icu_intime
     AND le.charttime <  b.icu_intime + INTERVAL '24 hours'
)
SELECT stay_id,
       MAX(CASE WHEN itemid = 50863 AND rn = 1 THEN valuenum END) AS alp,
       MAX(CASE WHEN itemid = 50863 AND rn = 1 THEN charttime END) AS alp_time,
       MAX(CASE WHEN itemid = 50862 AND rn = 1 THEN valuenum END) AS albumin,
       MAX(CASE WHEN itemid = 50862 AND rn = 1 THEN charttime END) AS albumin_time,
       MAX(CASE WHEN itemid = 50861 AND rn = 1 THEN valuenum END) AS alt,
       MAX(CASE WHEN itemid = 50878 AND rn = 1 THEN valuenum END) AS ast,
       MAX(CASE WHEN itemid = 50885 AND rn = 1 THEN valuenum END) AS bilirubin_total,
       MAX(CASE WHEN itemid = 50883 AND rn = 1 THEN valuenum END) AS bilirubin_direct,
       MAX(CASE WHEN itemid = 50912 AND rn = 1 THEN valuenum END) AS creatinine,
       MAX(CASE WHEN itemid = 50983 AND rn = 1 THEN valuenum END) AS sodium,
       MAX(CASE WHEN itemid = 50813 AND rn = 1 THEN valuenum END) AS lactate,
       MAX(CASE WHEN itemid = 51274 AND rn = 1 THEN valuenum END) AS pt,
       MAX(CASE WHEN itemid = 51237 AND rn = 1 THEN valuenum END) AS inr,
       MAX(CASE WHEN itemid = 51222 AND rn = 1 THEN valuenum END) AS hemoglobin,
       MAX(CASE WHEN itemid = 51265 AND rn = 1 THEN valuenum END) AS platelets,
       MAX(CASE WHEN itemid = 51301 AND rn = 1 THEN valuenum END) AS wbc,
       MAX(CASE WHEN itemid = 51279 AND rn = 1 THEN valuenum END) AS rbc,
       MAX(CASE WHEN itemid = 52075 AND rn = 1 THEN valuenum END) AS anc,
       MAX(CASE WHEN itemid = 51133 AND rn = 1 THEN valuenum END) AS alc,
       MAX(CASE WHEN itemid = 50852 AND rn = 1 THEN valuenum END) AS hba1c,
       MAX(CASE WHEN itemid = 50931 AND rn = 1 THEN valuenum END) AS glucose,
       MAX(CASE WHEN itemid = 50820 AND rn = 1 THEN valuenum END) AS ph,
       MAX(CASE WHEN itemid = 50818 AND rn = 1 THEN valuenum END) AS pco2,
       MAX(CASE WHEN itemid = 50821 AND rn = 1 THEN valuenum END) AS po2,
       MAX(CASE WHEN itemid = 50882 AND rn = 1 THEN valuenum END) AS bicarbonate,
       MAX(CASE WHEN itemid = 51006 AND rn = 1 THEN valuenum END) AS bun,
       MAX(CASE WHEN itemid = 50963 AND rn = 1 THEN valuenum END) AS ntprobnp,
       MAX(CASE WHEN itemid = 50907 AND rn = 1 THEN valuenum END) AS cholesterol_total,
       MAX(CASE WHEN itemid = 50904 AND rn = 1 THEN valuenum END) AS cholesterol_hdl,
       MAX(CASE WHEN itemid = 50905 AND rn = 1 THEN valuenum END) AS cholesterol_ldl,
       MAX(CASE WHEN itemid = 51000 AND rn = 1 THEN valuenum END) AS triglycerides
FROM ranked
GROUP BY stay_id;

/* ============================================================================
   STEP 5. EXCLUSION 1 - no valid ALP / albumin measurement WITHIN 24 h
   ----------------------------------------------------------------------------
   1,141 of the 3,068 patients with cirrhosis and sepsis had no ALP or no
   albumin measured in the first 24 h of the ICU stay and were excluded.
   ============================================================================ */
DROP TABLE IF EXISTS t_with_apar;
CREATE TEMP TABLE t_with_apar AS
SELECT s.*,
       l.alp, l.albumin, l.alp_time, l.albumin_time,
       l.alt, l.ast, l.bilirubin_total, l.bilirubin_direct,
       l.creatinine, l.sodium, l.lactate, l.pt, l.inr,
       l.hemoglobin, l.platelets, l.wbc, l.rbc, l.anc, l.alc,
       l.hba1c, l.glucose, l.ph, l.pco2, l.po2, l.bicarbonate,
       l.bun, l.ntprobnp,
       l.cholesterol_total, l.cholesterol_hdl, l.cholesterol_ldl, l.triglycerides,
       -- ALP (IU/L) / albumin (g/L); `albumin` below is the stored g/dL value
       (l.alp / NULLIF(l.albumin * 10.0, 0)) AS apar,
       (l.albumin * 10.0)                    AS albumin_gl
FROM t_sepsis s
JOIN t_lab l ON l.stay_id = s.stay_id
WHERE l.alp IS NOT NULL
  AND l.albumin IS NOT NULL
  AND l.albumin > 0;

-- expected: 1,927 (1,141 patients excluded)
SELECT COUNT(*) AS n_with_apar FROM t_with_apar;


/* ============================================================================
   STEP 6. EXCLUSION 2 - death within 24 h of ICU admission
   ----------------------------------------------------------------------------
   Patients who died within 24 h of ICU admission were excluded because the
   exposure cannot be interpreted as a baseline prognostic marker in patients
   who die before any clinical course can unfold.  Death is taken from
   admissions.deathtime (in-hospital) or patients.dod (death registry),
   whichever is recorded.

   NOTE ON ORDER: this exclusion is applied AFTER the laboratory exclusion and
   the 3-SD exclusion is NOT applied at all in the primary analysis (see the
   header).  Retaining these 62 patients is reported as sensitivity analysis
   11 (Appendix D).
   ============================================================================ */
DROP TABLE IF EXISTS t_alive24;
CREATE TEMP TABLE t_alive24 AS
SELECT w.*
FROM t_with_apar w
JOIN mimiciv_hosp.admissions ad ON ad.hadm_id = w.hadm_id
JOIN mimiciv_hosp.patients  p  ON p.subject_id = w.subject_id
WHERE COALESCE(ad.deathtime, p.dod::timestamp) IS NULL
   OR COALESCE(ad.deathtime, p.dod::timestamp) > w.icu_intime + INTERVAL '24 hours';

-- expected: 1,865 (62 patients excluded)
SELECT COUNT(*) AS n_after_24h_exclusion FROM t_alive24;


/* ============================================================================
   STEP 7. OUTCOMES - 30- and 90-day all-cause mortality
   ----------------------------------------------------------------------------
   Death dates in MIMIC-IV v3.1 (`patients.dod`) combine hospital records with
   state death-registry records, so deaths occurring after discharge are
   captured.  Patients are administratively censored at 30 and 90 days after
   ICU admission; because `dod` is a date (not a timestamp), follow-up time is
   counted in whole days.
   ============================================================================ */
DROP TABLE IF EXISTS t_outcome;
CREATE TEMP TABLE t_outcome AS
SELECT n.*,
       p.dod,
       (p.dod - n.icu_intime::date) AS days_to_death,
       CASE WHEN p.dod IS NOT NULL AND (p.dod - n.icu_intime::date) <= 30
            THEN 1 ELSE 0 END AS death_30d,
       CASE WHEN p.dod IS NOT NULL AND (p.dod - n.icu_intime::date) <= 90
            THEN 1 ELSE 0 END AS death_90d,
       CASE WHEN p.dod IS NULL THEN 30
            WHEN (p.dod - n.icu_intime::date) > 30 THEN 30
            ELSE (p.dod - n.icu_intime::date) END AS survival_days_30,
       CASE WHEN p.dod IS NULL THEN 90
            WHEN (p.dod - n.icu_intime::date) > 90 THEN 90
            ELSE (p.dod - n.icu_intime::date) END AS survival_days_90
FROM t_alive24 n
JOIN mimiciv_hosp.patients p ON p.subject_id = n.subject_id;


/* ============================================================================
   STEP 8. SEVERITY SCORES, AETIOLOGY, AKI, ACLF
   ----------------------------------------------------------------------------
   * SOFA, APS III, OASIS, SAPS II, GCS: official first-day derived
     values; their window runs from 6 h before to 24 h after ICU
     admission, which begins 6 h earlier than the exposure window.
   * MELD / MELD-Na: computed from the 24-h creatinine, bilirubin, sodium and
     the laboratory-reported INR, using the FULL set of 2016 OPTN
     bounding rules (OPTN Policy 9, effective January 2016):

         Cr_c     = min( max(creatinine, 1.0), 4.0 )   -- floor 1.0, ceiling 4.0
         Bili_c   = max(bilirubin, 1.0)                -- floor 1.0
         INR_c    = max(INR, 1.0)                      -- floor 1.0

     INR is taken directly from the laboratory report (MIMIC-IV item 51237),
     NOT approximated as PT / 11.0 s.  The approximation would require both
     an international sensitivity index of 1 and a normal mean prothrombin
     time of exactly 11.0 s; in this cohort the reported INR implies a normal
     mean PT with a median of 10.889 s, the two definitions correlate at
     r = 0.9722 and agree exactly in only 8.4% of patients.  The difference is
     not a timing artefact: PT and INR are drawn from the same charttime in
     1,831 of 1,833 patients (99.89%).  Using the reported INR also keeps
     MELD-Na consistent with the ACLF coagulation criterion below, which uses
     the same item.  Prothrombin time is still extracted (item 51274) and
     reported in Table 1; the PT / 11.0 s variant is retained as a sensitivity
     analysis in Supplementary Table S5.
         MELD(i)  = max( 0.957*ln(Cr_c) + 0.378*ln(Bili_c)
                         + 1.120*ln(INR_c) + 0.643 , 0 ) * 10
         Na_c     = min(137, max(125, sodium))    -- OPTN sodium window 125-137
                  = LEAST(GREATEST(sodium, 125.0), 137.0)
         sodium term (applied ONLY when MELD(i) > 11):
                    1.32*(137-Na_c) - 0.033*MELD(i)*(137-Na_c)
         MELD-Na  = min( max( MELD(i) + sodium_term , 0 ), 40 )

     Both sodium endpoints are inclusive: Na = 120 -> 125, Na = 125 -> 125,
     Na = 130 -> 130, Na = 137 -> 137, Na = 140 -> 137.  The clamp is written
     as a continuous LEAST/GREATEST pair, so no observed value can land in an
     undefined interval.  The identical rules and the same test values are
     asserted at import time by _self_test_meld_bounds() in apar_common.py
     (Supplementary File 3), which is the single implementation used by every
     analysis script.  MELD(i) is kept continuous rather than rounded to an
     integer because it is used as a continuous covariate.
     Rows with a missing component yield a NULL MELD-Na and are handled by the
     complete-case rule of the analysis.  In this cohort 32 patients lack the
     INR (1.7%; the missingness pattern of INR and PT is identical) and a
     further 2 lack bilirubin, so neither score can be computed for 34 patients
     (1.8%) - this is the origin of the 1,831 complete cases.
   * Aetiology: alcohol-related cirrhosis = ICD-9 5712 or ICD-10 K7030/K7031/
     K7032 in any diagnosis position.
   * AKI: ICD-9 584.5-584.9 or ICD-10 N17.0-N17.9 in any diagnosis position.
     Codes are compared after the decimal point is removed and matched as
     prefixes ('584%', 'N17%'), so N17.8 and N17.9 are included; see
     Supplementary Table S10.
   * ACLF: >= 1 organ failure inside the first 24 h, using the six components
     listed below (simplified CLIF-C Organ Failure definition):

         liver         bilirubin >= 12.0 mg/dL (24-h maximum)
         kidney        creatinine >= 2.0 mg/dL (24-h maximum) OR renal
                       replacement therapy within 24 h
         coagulation   INR >= 2.5 (24-h maximum; item 51237, the same
                       laboratory-reported INR used by MELD-Na)
         circulation   any vasopressor within 24 h
         brain         GCS <= 8 (minimum over the derived first-day window,
                       taken from the official table
                       mimiciv_derived.first_day_gcs, which aggregates from
                       6 h before to 24 h after ICU admission and is thus a
                       little wider than the 0-24 h window used by the other
                       five components; see Supplementary Table S10).  The
                       whole-stay minimum GCS is also
                       extracted, as the column gcs_whole, but it is used by
                       no analysis other than Table S5 analysis 15, which
                       deliberately re-imposes the wider window to quantify
                       the effect of this choice.
         respiration   PaO2/FiO2 <= 200, using the 24-h minimum PaO2 and the
                       24-h maximum FiO2 (worst-case pairing)

     FiO2 ITEM AND UNIT (corrected 2026-10-24):  FiO2 is taken from
         chartevents item 223835 (Inspired O2 Fraction) only.  An earlier
         version of this script also listed items 220210, 3420 and 3422;
         220210 is "Respiratory Rate", not FiO2 (mimiciv_icu.d_items), and
         3420 / 3422 are CareVue items that carry no rows in MIMIC-IV, so the
         earlier MAX() pooled respiratory rates (53,934 rows, 1,862 stays)
         with true FiO2 values (7,412 rows, 1,142 stays).  A recorded value
         is retained only when it is physiologically interpretable:
                 0 < v <= 1.5     recorded as a fraction, used as it stands
                 21 <= v <= 100   recorded as a percentage, divided by 100
                 otherwise        discarded (v <= 0, 1.5 < v < 21, v > 100)
         Values above 1.5 are percentages and are divided by 100 before the
         ratio is formed; values at or below 1.5 are fractions and are not.
         After the filter 1,141 of 1,865 stays (61.2%) have a usable FiO2 and
         the respiratory criterion is met by 927 (49.7%) of the cohort.
   ============================================================================ */
/* MELD(i) and MELD-Na are built in their own temp table so that the
   "sodium term only when MELD(i) > 11" rule can reference MELD(i) without
   repeating the whole logarithm expression three times.

   PRECISION NOTE: both scores are kept continuous and are NOT rounded, so that
   the columns written to apar_cohort_1865 are identical, row by row, to the
   values that apar_common.meld_from() / meld_na_from() recompute from the same
   raw laboratory values in Supplementary File 3 (maximum absolute difference
   over the 1,865 patients: 0.000 for MELD(i)).  MELD(i) is exported as `meld`
   so that the two liver-specific rows of Table 1 can be reproduced from the
   final table alone. */
DROP TABLE IF EXISTS t_meld;
CREATE TEMP TABLE t_meld AS
WITH b AS (
  SELECT o.stay_id,
         o.sodium,
         CASE WHEN o.creatinine IS NULL OR o.creatinine <= 0
                OR o.bilirubin_total IS NULL OR o.bilirubin_total <= 0
                OR o.inr IS NULL OR o.inr <= 0
              THEN NULL
              ELSE GREATEST(
                     0.957 * LN(LEAST(GREATEST(o.creatinine, 1.0), 4.0))
                   + 0.378 * LN(GREATEST(o.bilirubin_total, 1.0))
                   + 1.120 * LN(GREATEST(o.inr, 1.0))
                   + 0.643, 0) * 10
         END AS meld_i
  FROM t_outcome o
)
SELECT b.stay_id,
       b.meld_i AS meld,
       CASE
         WHEN b.meld_i IS NULL OR b.sodium IS NULL THEN NULL
         ELSE LEAST(GREATEST(
                  b.meld_i
                + CASE WHEN b.meld_i > 11.0
                       THEN 1.32 * (137.0 - LEAST(GREATEST(b.sodium, 125.0), 137.0))
                          - 0.033 * b.meld_i
                            * (137.0 - LEAST(GREATEST(b.sodium, 125.0), 137.0))
                       ELSE 0.0 END, 0), 40)
       END AS meld_na
FROM b;

DROP TABLE IF EXISTS t_scores;
CREATE TEMP TABLE t_scores AS
SELECT o.*,
       fds.sofa  AS sofa,
       ap.apsiii AS apsiii,
       oa.oasis  AS oasis,
       sp.sapsii AS sapsii,
       gc.gcs_min AS gcs,
       -- the minimum over the whole stay is kept alongside it, but is used
       -- only by Table S5 analysis 15 (see the note at the ACLF block above)
       gcw.gcs_min AS gcs_whole,
       ch.charlson_comorbidity_index AS charlson,
       -- MELD (MELD(i)) and MELD-Na from the 24-h values, full 2016 OPTN rules
       -- (see the note above)
       tm.meld     AS meld,
       tm.meld_na  AS meld_na
FROM t_outcome o
LEFT JOIN t_meld tm ON tm.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.first_day_sofa fds ON fds.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.apsiii  ap ON ap.stay_id  = o.stay_id
LEFT JOIN mimiciv_derived.oasis   oa ON oa.stay_id  = o.stay_id
LEFT JOIN mimiciv_derived.sapsii  sp ON sp.stay_id  = o.stay_id
LEFT JOIN mimiciv_derived.first_day_gcs gc ON gc.stay_id = o.stay_id
LEFT JOIN (SELECT stay_id, MIN(gcs) AS gcs_min FROM mimiciv_derived.gcs
           GROUP BY stay_id) gcw ON gcw.stay_id = o.stay_id
LEFT JOIN mimiciv_derived.charlson ch ON ch.hadm_id = o.hadm_id;

-- aetiology (alcohol-related cirrhosis) and AKI
DROP TABLE IF EXISTS t_dx;
CREATE TEMP TABLE t_dx AS
SELECT o.stay_id,
       MAX(CASE WHEN (d.icd_version = 9  AND REPLACE(d.icd_code, '.', '') = '5712')
                  OR (d.icd_version = 10 AND UPPER(REPLACE(d.icd_code, '.', ''))
                      IN ('K7030','K7031','K7032'))
                THEN 1 ELSE 0 END) AS alcohol_flag,
       MAX(CASE WHEN (d.icd_version = 9  AND REPLACE(d.icd_code, '.', '') LIKE '584%')
                  OR (d.icd_version = 10 AND UPPER(REPLACE(d.icd_code, '.', '')) LIKE 'N17%')
                THEN 1 ELSE 0 END) AS aki,

       -- biliary obstruction / cholangitis / gallstone disease.
       -- The code list is identical to Supplementary File 3,
       -- sql/03_sensitivity_cohorts.sql BLOCK 4, so that this flag reproduces
       -- the same 152 cohort patients that sensitivity analysis 6 starts from
       -- (Table S5 row 6): 150 of them have complete Model 3 covariates, so
       -- the analysis is fitted on 1,831 - 150 = 1,681.  ICD-10 codes are
       -- prefix matches: the bare codes
       -- 'K80' and 'K81' never occur in MIMIC-IV (only K8000-K8081 and
       -- K810 / K811 / K812 / K819), whereas 'K830' and 'K831' do occur and
       -- K830 additionally has K8301 / K8309.  Any diagnosis position.
       MAX(CASE WHEN (d.icd_version = 9  AND REPLACE(d.icd_code, '.', '') IN
                        ('5761','5762','5740','5741','5742','5743','5744','5745',
                         '5746','5747','5748','5749','5750'))
                  OR (d.icd_version = 10 AND UPPER(REPLACE(d.icd_code, '.', ''))
                        LIKE ANY (ARRAY['K830%','K831%','K80%','K81%']))
                THEN 1 ELSE 0 END) AS biliary_obs
FROM t_scores o
JOIN mimiciv_hosp.diagnoses_icd d
  ON d.subject_id = o.subject_id AND d.hadm_id = o.hadm_id
GROUP BY o.stay_id;

-- ACLF: organ failures inside the first 24 h
DROP TABLE IF EXISTS t_aclf;
CREATE TEMP TABLE t_aclf AS
WITH w AS (SELECT stay_id, hadm_id, icu_intime,
                  icu_intime + INTERVAL '24 hours' AS t24 FROM t_scores),
agg AS (
  SELECT w.stay_id,
         MAX(CASE WHEN l.itemid = 50885 THEN l.valuenum END) AS bili_max,
         MAX(CASE WHEN l.itemid = 50912 THEN l.valuenum END) AS cr_max,
         MAX(CASE WHEN l.itemid = 51237 THEN l.valuenum END) AS inr_max,
         MIN(CASE WHEN l.itemid = 50821 THEN l.valuenum END) AS pao2_min
  FROM w
  JOIN mimiciv_hosp.labevents l
    ON l.hadm_id = w.hadm_id
   AND l.itemid IN (50885, 50912, 51237, 50821)
   AND l.valuenum IS NOT NULL
   AND l.charttime >= w.icu_intime AND l.charttime < w.t24
  GROUP BY w.stay_id),
rrt  AS (SELECT p.stay_id, 1 AS rrt_flag FROM mimiciv_icu.procedureevents p
         JOIN w ON w.stay_id = p.stay_id
         WHERE p.starttime >= w.icu_intime AND p.starttime < w.t24
           AND p.itemid IN (225802,225803,225809,225955,225441,225805,225436)
         GROUP BY p.stay_id),
vaso AS (SELECT i.stay_id, 1 AS vaso_flag FROM mimiciv_icu.inputevents i
         JOIN w ON w.stay_id = i.stay_id
         WHERE i.starttime >= w.icu_intime AND i.starttime < w.t24
           AND i.itemid IN (221906,221289,229617,221749,229632,229789,222315,221662)
         GROUP BY i.stay_id),
fio2 AS (SELECT c.stay_id, MAX(c.valuenum) AS fio2_max FROM mimiciv_icu.chartevents c
         JOIN w ON w.stay_id = c.stay_id
         WHERE c.charttime >= w.icu_intime AND c.charttime < w.t24
           AND c.itemid = 223835
           AND c.valuenum IS NOT NULL
           AND ( (c.valuenum > 0   AND c.valuenum <= 1.5)
              OR (c.valuenum >= 21 AND c.valuenum <= 100) )
         GROUP BY c.stay_id)
SELECT w.stay_id,
       COALESCE(agg.bili_max >= 12.0, FALSE)::int                       AS of_liver,
       (COALESCE(agg.cr_max >= 2.0, FALSE) OR COALESCE(rrt.rrt_flag,0)=1)::int AS of_kidney,
       COALESCE(agg.inr_max >= 2.5, FALSE)::int                         AS of_coagulation,
       COALESCE(vaso.vaso_flag,0)                                       AS of_circulation,
       COALESCE(fio2.fio2_max, 0)                                       AS fio2_max,
       COALESCE(agg.pao2_min, NULL)                                     AS pao2_min
FROM w
LEFT JOIN agg  USING (stay_id)
LEFT JOIN rrt  USING (stay_id)
LEFT JOIN vaso USING (stay_id)
LEFT JOIN fio2 USING (stay_id);


/* ============================================================================
   STEP 9. FINAL ANALYTIC TABLE
   ============================================================================ */
DROP TABLE IF EXISTS apar_cohort_1865;
CREATE TABLE apar_cohort_1865 AS
SELECT s.subject_id,
       s.hadm_id,
       s.stay_id,
       s.icu_intime,
       s.icu_outtime,
       s.los_icu,
       EXTRACT(EPOCH FROM (s.sepsis_onset_time - s.icu_intime)) / 3600.0
           AS hours_icu_to_sepsis_onset,
       s.suspected_infection_time,
       s.sepsis_onset_time,
       s.sofa_at_onset,
       s.respiration, s.coagulation, s.liver,
       s.cardiovascular, s.cns, s.renal,
       s.alp, s.albumin, s.albumin_gl, s.alp_time, s.albumin_time,
       EXTRACT(EPOCH FROM (s.alp_time      - s.icu_intime)) / 3600.0 AS alp_h_from_icu,
       EXTRACT(EPOCH FROM (s.albumin_time  - s.icu_intime)) / 3600.0 AS alb_h_from_icu,
       ROUND(s.apar::numeric, 4) AS apar,
       s.gender,
       s.age,
       s.sofa, s.apsiii, s.oasis, s.sapsii, s.gcs, s.gcs_whole, s.charlson,
       s.creatinine, s.bilirubin_total, s.sodium, s.pt, s.inr, s.lactate,
       s.hemoglobin, s.platelets, s.wbc, s.alt, s.ast,
       s.meld, s.meld_na,
       COALESCE(x.alcohol_flag, 0) AS alcohol_flag,
       COALESCE(x.aki, 0)          AS aki,
       COALESCE(x.biliary_obs, 0)  AS biliary_obs,
       ( COALESCE(a.of_liver, 0) + COALESCE(a.of_kidney, 0)
       + COALESCE(a.of_coagulation, 0) + COALESCE(a.of_circulation, 0)
       + (s.gcs IS NOT NULL AND s.gcs <= 8)::int
       + (a.pao2_min IS NOT NULL AND a.fio2_max > 0
          AND a.pao2_min / (CASE WHEN a.fio2_max > 1.5 THEN a.fio2_max / 100.0
                                 ELSE a.fio2_max END) <= 200)::int ) AS of_count,
       s.dod,
       s.days_to_death,
       s.death_30d,
       s.death_90d,
       s.survival_days_30,
       s.survival_days_90
FROM t_scores s
LEFT JOIN t_dx   x ON x.stay_id = s.stay_id
LEFT JOIN t_aclf a ON a.stay_id = s.stay_id;

ALTER TABLE apar_cohort_1865 ADD COLUMN aclf_binary int;
UPDATE apar_cohort_1865 SET aclf_binary = (of_count >= 1)::int;


/* ============================================================================
   STEP 10. VERIFICATION - these numbers must match the manuscript
   ============================================================================ */
SELECT COUNT(*)                                AS n_cohort     FROM apar_cohort_1865;  -- 1,865
SELECT COUNT(*) FILTER (WHERE death_30d = 1)   AS n_death_30d  FROM apar_cohort_1865;  --   599
SELECT COUNT(*) FILTER (WHERE death_90d = 1)   AS n_death_90d  FROM apar_cohort_1865;  --   755
SELECT ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY apar)::numeric, 2) AS apar_q1,
       ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY apar)::numeric, 2) AS apar_median,
       ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY apar)::numeric, 2) AS apar_q3,
       ROUND(AVG(apar)::numeric, 2)            AS apar_mean,
       ROUND(STDDEV_SAMP(apar)::numeric, 2)    AS apar_sd
FROM apar_cohort_1865;                       -- 2.43 / 3.46 / 5.39 ; mean 4.77 ; SD 4.71
WITH s AS (SELECT AVG(apar) AS m, STDDEV_SAMP(apar) AS sd FROM apar_cohort_1865)
SELECT ROUND((s.m + 3 * s.sd)::numeric, 2) AS apar_3sd_threshold,
       (SELECT COUNT(*) FROM apar_cohort_1865 c WHERE c.apar > s.m + 3 * s.sd) AS n_above_3sd
FROM s;                                      -- 18.89 ; 36 (RETAINED in the primary analysis)

-- Albumin is reported in g/L in the manuscript; the stored value is g/dL.
SELECT ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY albumin_gl)::numeric, 1)
           AS albumin_g_l_median
FROM apar_cohort_1865;                       -- 29.0

-- Biliary obstruction / cholangitis / gallstone disease
-- (identical code list to Supplementary File 3, BLOCK 4; sensitivity 6)
SELECT COUNT(*) FILTER (WHERE biliary_obs = 1) AS n_biliary
FROM apar_cohort_1865;                   -- 152 (cohort); 150 of them are
                                         -- complete cases, so analysis 6
                                         -- is fitted on 1,831 - 150 = 1,681

-- MELD and MELD-Na reproduce the two liver-specific rows of Table 1.  Both
-- come from the 24-h creatinine, bilirubin and sodium and the laboratory-
-- reported INR; 34 of the 1,865 patients have no computable score (32 lack the
-- INR, 2 lack bilirubin), which is the origin of the 1,831 complete cases.
SELECT COUNT(meld)                          AS n_meld_avail FROM apar_cohort_1865;  -- 1,831
SELECT ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY meld)::numeric, 1) AS meld_q1,
       ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY meld)::numeric, 1) AS meld_median,
       ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY meld)::numeric, 1) AS meld_q3
FROM apar_cohort_1865;                       -- 15.5 / 21.5 / 28.8
SELECT ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY meld_na)::numeric, 1) AS meld_na_q1,
       ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY meld_na)::numeric, 1) AS meld_na_median,
       ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY meld_na)::numeric, 1) AS meld_na_q3
FROM apar_cohort_1865;                       -- 16.3 / 23.0 / 30.5

-- ACLF (EASL-CLIF) and its brain-failure component.  The brain criterion
-- uses the minimum GCS of the official first-day derived table, whose
-- window runs from 6 h before to 24 h after ICU admission; the whole-stay
-- minimum is kept as gcs_whole and is used only
-- by Table S5 analysis 15.  Expected counts are the manuscript values,
-- re-verified against mimiciv_hosp.apar_cohort_1865 on 2026-10-25; they
-- supersede those carried by the pre-FiO2-correction version of this
-- script (1,549 patients with ACLF, grades 488 / 490 / 571).  Grades sum
-- to the total: 517 + 462 + 516 = 1,495.
SELECT COUNT(*) FILTER (WHERE aclf_binary = 1)      AS n_aclf,             -- 1,495
       COUNT(*) FILTER (WHERE gcs IS NOT NULL AND gcs <= 8)
                                                    AS n_brain_24h,        --   157
       COUNT(*) FILTER (WHERE gcs_whole IS NOT NULL AND gcs_whole <= 8)
                                                    AS n_brain_wholestay,  --   413
       COUNT(*) FILTER (WHERE of_count = 1)         AS n_aclf_grade1,      --   517
       COUNT(*) FILTER (WHERE of_count = 2)         AS n_aclf_grade2,      --   462
       COUNT(*) FILTER (WHERE of_count >= 3)        AS n_aclf_grade3       --   516
FROM apar_cohort_1865;

-- Export (optional)
-- \copy apar_cohort_1865 TO 'apar_cohort_1865.csv' CSV HEADER;


/* ============================================================================
   APPENDIX A. SUSPECTED INFECTION - the window, as implemented in
               mimiciv_derived.sepsis3
   ----------------------------------------------------------------------------
   A patient is considered to have a suspected infection if BOTH of the
   following occur within a single hospital admission:

     (a) a new systemic antibacterial agent is started
         (mimiciv_hosp.prescriptions, identified by generic drug name), AND
     (b) a clinical microbiological specimen is obtained
         (mimiciv_hosp.microbiologyevents)

   WINDOW: the two events above must additionally fall close together in
   time.  EITHER of the following is sufficient - they are alternatives,
   not requirements that must hold simultaneously:

     (i)  the antibacterial agent is started no more than 72 h AFTER the
          specimen is obtained, OR
     (ii) the specimen is obtained no more than 24 h AFTER the
          antibacterial agent is started

   The union of (i) and (ii) is the single interval

       culture_time BETWEEN antibiotic_time - 72 h
                        AND antibiotic_time + 24 h

   NOTE ON THE -72/+24 BOUNDS: the specimen may therefore be obtained up to
   72 h BEFORE, or up to 24 h AFTER, the first antibacterial dose.  The SQL
   above reproduces the two criteria implemented in the official source
   distributed with the MIMIC Code Repository as
       mimic-iv/concepts/sepsis/suspicion_of_infection.sql
   (CTEs me_then_ab and ab_then_me).  Where a specimen carries no timestamp,
   that source falls back to a date-level match over 3 days, which is why a
   small residual number of pairs lie marginally outside -72 h.

   ANCHOR (suspected_infection_time): the earlier of the two times
       suspected_infection_time = LEAST(antibiotic_time, culture_time)
       (if only one of the two is present, that time is used)

   SEPSIS-3 ONSET (t0): the first time at which SOFA >= 2 within the window
       suspected_infection_time - 48 h  ..  suspected_infection_time + 24 h,
       with the additional requirement that the SOFA increase be attributable
       to infection (operationalised as SOFA >= 2 in that window).
   ============================================================================ */

/* ============================================================================
   APPENDIX B. SOFA - calculation method and baseline definition
   ----------------------------------------------------------------------------
   The SOFA score that operationalises the Sepsis-3 criterion here is the
   total SOFA at sepsis onset (t0), taken from
   `mimiciv_derived.sepsis3.sofa_score`; its six sub-scores are stored in the
   same table (respiration ... renal). It is NOT the covariate used in the
   models - the two SOFA columns are contrasted under "SOFA USED IN THE
   MANUSCRIPT" at the end of this appendix.  The sub-scores are computed by
   the official MIMIC-IV derived table `mimiciv_derived.sofa` (Guan et al.),
   which scores the six organ systems hour by hour over the ICU stay.  The
   components and the inputs each of them uses are:

     Respiration             PaO2/FiO2 ratio, scored through two separate
                             channels according to invasive ventilation
                             status
     Coagulation             platelet count
     Liver                   total bilirubin
     Cardiovascular          mean arterial pressure; infusion rates of
                             dopamine, dobutamine, epinephrine and
                             norepinephrine
     Central nervous system  Glasgow Coma Scale
     Renal                   creatinine; urine output

   The inputs are taken, as in the official code, from
       mimiciv_icu.chartevents    GCS, mean arterial pressure, FiO2,
                                  ventilation status
       mimiciv_hosp.labevents     creatinine, bilirubin, PaO2, platelets
       mimiciv_icu.inputevents    vasopressor infusions
       mimiciv_icu.outputevents   urine output
   Each hourly sub-score is subsequently replaced by its maximum over the
   preceding 24 h of the ICU stay (see below), and the six maxima are summed
   to give the total.

   The per-grade cut-points are NOT reproduced here.  They are defined by the
   official code distributed with the MIMIC Code Repository as
       mimic-iv/concepts/score/sofa.sql
   and are applied verbatim by `mimiciv_derived.sofa`; that source - not any
   secondary summary of it - is authoritative for the exact thresholds.  Two
   of the six components in particular resist a compact tabulation: the
   respiration component is scored separately for ventilated and
   non-ventilated measurements, and the renal component combines creatinine
   with urine output instead of following a single ladder of creatinine
   bands.  Both are implemented exactly as in the official code.

   BASELINE SOFA AND HOW IT ENTERS THE SEPSIS-3 CRITERION
     Sepsis-3 asks for an ACUTE INCREASE of >= 2 SOFA points attributable to
     infection, which strictly presupposes knowledge of pre-morbid organ
     function. Neither `mimiciv_derived.sofa` nor `mimiciv_derived.sepsis3`
     performs such a baseline correction: both score organ function de novo.
     The official `concepts/sepsis/sepsis3.sql` says so explicitly -

        "-- extract rows with SOFA >= 2
         -- implicitly this assumes baseline SOFA was 0 before ICU admission."

        "Implicitly, the baseline SOFA score is assumed to be zero, as we do
         not know if the patient has preexisting (acute or chronic) organ
         dysfunction before the onset of infection."

     - so the required >= 2-point increase is OPERATIONALISED as an absolute
     total SOFA of >= 2, with baseline SOFA assumed to be 0. Every primary
     analysis in this study uses that definition, exactly as the manuscript
     states in Methods 2.2. Verified locally: all 41,295 rows flagged
     `sepsis3 = TRUE` in mimiciv_derived.sepsis3 have sofa_score >= 2
     (minimum = 2), i.e. the derived table returns an absolute value of >= 2
     and not a change from baseline.

     `mimiciv_derived.sofa` does refer to a 24-hour window, but for a
     different purpose and against a different anchor: each hourly sub-score
     is replaced by its maximum over the preceding 24 h of the ICU stay,
     MAX(component) OVER (PARTITION BY stay_id ORDER BY hr
     ROWS BETWEEN 23 PRECEDING AND 0 FOLLOWING), and a component missing
     throughout is coalesced to 0. That rolling window is a carry-forward
     rule applied inside the ICU stay; it is NOT a pre-admission baseline
     and must not be read as one.

   SOFA USED IN THE MANUSCRIPT
     * Covariate (Methods 2.4, every model): `fds.sofa`, i.e.
       `mimiciv_derived.first_day_sofa.sofa` - the first-day score, built
       from the worst value of each of the six components over a window
       running from 6 h before to 24 h after ICU admission; the per-grade
       cut-points are those of mimic-iv/concepts/firstday/first_day_sofa.sql.
       This is the column carried into apar_cohort_1865 as `sofa`.
     * Sepsis-3 definition (Methods 2.2): `sofa_score` at `sofa_time` (t0) from
       mimiciv_derived.sepsis3 - the SOFA at sepsis onset.  The two are scored
       by the same six-component rule but over different windows, so they are
       reported separately rather than used interchangeably; in this cohort
       they differ in 1,645 of the 1,865 patients (median 9 vs 5).
     * Sensitivity analysis 7 (Supplementary Table S5, row 7; numbered 7 in
       this file and in Supplementary File 3, script 07_sensitivity.py):
       sepsis re-defined by a genuinely baseline-corrected RISE.
           baseline SOFA = MIN(sofa_24hours) over the 48 h preceding
                           suspected_infection_time (anchor = suspected
                           infection, NOT ICU admission)
           SOFA at onset = the sofa_24hours row closest to
                           sepsis3.sofa_time (within +/- 30 min)
           dSOFA         = SOFA at onset - baseline SOFA
       A patient is excluded only when a baseline IS available and the rise
       is < 2; the 1,102 patients with no computable baseline (suspected
       infection preceded ICU admission, so the criterion cannot be shown to
       have failed rather than is unmeasured) are RETAINED. In this cohort
       382 patients reached dSOFA >= 2, 381 were excluded, leaving 1,484
       patients, of whom 1,459 had complete Model 3 covariates. The PRIMARY
       analyses do not use this correction - see BASELINE SOFA above.
   ============================================================================ */

/* ============================================================================
   APPENDIX C. SAMPLING-TIME DISTRIBUTION (Supplementary Table S2)
   ----------------------------------------------------------------------------
   Because ALP and albumin may be drawn at different times, the time elapsed
   between ICU admission and each index measurement, and the delay between the
   two draws, are reported in Supplementary Table S2.  The query below
   reproduces that table from the analytic cohort.

   Reported values (N = 1,865):
     ALP      : median 1.80 h after ICU admission (IQR 0.61 - 6.18)
     Albumin  : median 2.38 h after ICU admission (IQR 0.69 - 7.97)
     Same tube / same draw (identical charttime): 1,645 patients (88.2%)
     Albumin drawn later than ALP: 220 patients (11.8%), median delay 9.76 h
   ============================================================================ */
SELECT COUNT(*) FILTER (WHERE alp_time = albumin_time)                     AS n_same_draw,
       ROUND(100.0 * COUNT(*) FILTER (WHERE alp_time = albumin_time)
             / COUNT(*), 1)                                                AS pct_same_draw,
       ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY alp_h_from_icu)::numeric, 2) AS alp_median_h,
       ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY alb_h_from_icu)::numeric, 2) AS alb_median_h
FROM apar_cohort_1865;


/* ============================================================================
   APPENDIX D. SENSITIVITY COHORTS
   ----------------------------------------------------------------------------
   D1. APAR > mean + 3 SD excluded (sensitivity analysis 1)
       Threshold recomputed inside the cohort: mean + 3 SD = 18.89, n = 36
       excluded -> retained cohort N = 1,829
       (577 deaths at 30 days, 726 at 90 days).
       The model is fitted on the 1,795 complete cases
       (568 events at 30 days, 715 at 90 days); 34 patients have no
       computable MELD-Na.  The manuscript reports n = 1,795.

   D2. Deaths within 24 h RETAINED (sensitivity analysis 11)
       Omit STEP 6 (t_alive24 = t_with_apar) -> retained cohort N = 1,927
       (661 deaths at 30 days, 817 at 90 days).
       The model is fitted on the 1,891 complete cases
       (650 events at 30 days, 804 at 90 days); 36 patients have no
       computable MELD-Na.  The manuscript reports n = 1,891.

   Both variants are produced by the accompanying Python package
   (Supplementary File 3, script 07_sensitivity.py).
   ============================================================================ */

/* ============================================================================
   END OF FILE
   ============================================================================ */
