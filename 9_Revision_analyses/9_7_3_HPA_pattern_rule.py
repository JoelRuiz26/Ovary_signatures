#!/usr/bin/env python
"""Applies, per antibody, the rule 'absent or low in normal stroma and positive in most tumors'.

Rule (per antibody): normal ovarian stroma with staining in {not detected, low} AND
tumors with staining other than 'not detected' in more than 50 % of the patients.
The candidate status of each gene is read from 6_Depmap_ovary/6_3_2_candidate_regulators.tsv.

Input:  $OVARY_SIGNATURES/9_Revision_analyses/hpa/hpa_v25_ovary_by_antibody.tsv,
        .../hpa_v25_ovary_patients.tsv (9_7_2_HPA_extract.py)
Output: $OVARY_SIGNATURES/9_Revision_analyses/hpa/pattern_rule_by_antibody.tsv
Usage:  python 9_7_3_HPA_pattern_rule.py [by_antibody.tsv patients.tsv output.tsv]
"""
import os, sys, csv, collections

BASE = os.path.expanduser(os.environ.get("OVARY_SIGNATURES", "~/Ovary_signatures"))  # project data root
HPA_DIR = os.path.join(BASE, "9_Revision_analyses/hpa")
if len(sys.argv) >= 4:
    ab_tsv, pat_tsv, out = sys.argv[1:4]
else:
    ab_tsv = os.path.join(HPA_DIR, "hpa_v25_ovary_by_antibody.tsv")
    pat_tsv = os.path.join(HPA_DIR, "hpa_v25_ovary_patients.tsv")
    out = os.path.join(HPA_DIR, "pattern_rule_by_antibody.tsv")

# candidate regulators of the pipeline (NA if the gene was not evaluated there)
cand = {}
cand_tsv = os.path.join(BASE, "6_Depmap_ovary/6_3_2_candidate_regulators.tsv")
if os.path.exists(cand_tsv):
    cand = {r["gene"]: r["candidate"] for r in csv.DictReader(open(cand_tsv), delimiter="\t")}

# counts in serous tumors only (each patient counts once)
ser = collections.defaultdict(dict)
for p in csv.DictReader(open(pat_tsv), delimiter="\t"):
    if p["tissue"] == "Ovarian cancer" and "serous" in p["snomed"]:
        ser[(p["gene"], p["antibody"])][p["patientId"]] = p["staining"].lower()

rows = []
for r in csv.DictReader(open(ab_tsv), delimiter="\t"):
    if r["tumor_n_patients"] == "NA":
        continue
    n, det, mh = int(r["tumor_n_patients"]), int(r["tumor_n_staining_detected"]), int(r["tumor_n_staining_med_high"])
    stroma = r["normal_stroma_staining"]
    c1 = "NA" if stroma == "NA" else str(stroma in ("not detected", "low"))
    c2 = det / n > 0.5
    s = ser[(r["gene"], r["antibody"])]
    sdet = sum(v != "not detected" for v in s.values())
    rows.append(dict(gene=r["gene"], antibody=r["antibody"], stroma=stroma, follicle=r["normal_follicle_staining"],
                     n=n, detected=det, medium_high=mh, intensity_nonneg=r["tumor_n_intensity_nonneg"],
                     stroma_ND_or_low=c1, tumor_majority_detected=c2,
                     meets_rule=("NA" if c1 == "NA" else str(c1 == "True" and c2)),
                     majority_medium_high=mh / n > 0.5, serous_n=len(s), serous_detected=sdet,
                     candidate=cand.get(r["gene"], "NA")))
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for f in rows:
    print("\t".join(str(v) for v in f.values()))
missing = sorted(g for g, c in cand.items() if c == "TRUE" and g not in {r["gene"] for r in rows})
if missing:
    print("WARNING: candidates without HPA tumor IHC in the table:", ", ".join(missing), file=sys.stderr)
