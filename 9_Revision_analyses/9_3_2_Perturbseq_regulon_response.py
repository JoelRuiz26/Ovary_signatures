#!/usr/bin/env python
"""Transcriptional response to the CRISPRi knockdown of each candidate regulator
(Replogle et al. 2022; normalized pseudo-bulk = z-score against control cells)
compared with its ARACNe-AP regulon from TCGA-OV.

Test: the positive targets of the regulon go down more than the genes outside the
regulon (one-sided Mann-Whitney), on the raw z-score and on the residual after
removing the shared response (mean of all perturbations).
Specificity: rank of the regulator's own regulon among the regulons of all the
candidates in that knockdown.

Input:  9_Revision_analyses/perturbseq/enriched_regulons.tsv (9_3_1_Export_regulons.R)
        <dataset>_normalized_bulk_01.h5ad in REPLOGLE_DIR
        (default $OVARY_SIGNATURES/data/replogle2022)
Output: 9_Revision_analyses/perturbseq/perturbseq_regulon_tests.tsv
Usage:  python 9_3_2_Perturbseq_regulon_response.py [dataset ...]
        (default: K562_essential K562_gwps rpe1)
"""
import os, sys, numpy as np, pandas as pd, anndata as ad
from scipy.stats import mannwhitneyu, false_discovery_control

BASE = os.path.expanduser(os.environ.get("OVARY_SIGNATURES", "~/Ovary_signatures"))  # project data root
DATA = os.environ.get("REPLOGLE_DIR", os.path.join(BASE, "data/replogle2022"))
OUT_DIR = os.path.join(BASE, "9_Revision_analyses/perturbseq")
datasets = sys.argv[1:] or ["K562_essential", "K562_gwps", "rpe1"]

reg = pd.read_csv(os.path.join(OUT_DIR, "enriched_regulons.tsv"), sep="\t")
reg = reg[reg.candidate]
drivers = list(reg.driver.unique())
pos_sets = {d: set(reg[(reg.driver == d) & (reg.tfmode > 0)].target) for d in drivers}
all_sets = {d: set(reg[reg.driver == d].target) for d in drivers}

def mw(x, y):
    return mannwhitneyu(x, y, alternative="less").pvalue if len(x) > 2 else np.nan

rows = []
for ds in datasets:
    a = ad.read_h5ad(os.path.join(DATA, f"{ds}_normalized_bulk_01.h5ad"))
    X = pd.DataFrame(np.asarray(a.X), index=a.obs_names, columns=a.var["gene_name"].values)
    X = X.loc[:, ~X.columns.duplicated()].replace([np.inf, -np.inf], np.nan)
    X = X.loc[:, X.notna().all(axis=0)].astype("float64")
    pert_gene = a.obs_names.str.split("_").str[1]
    generic = X.mean(axis=0); g0 = generic - generic.mean()
    for d in drivers:
        idx = np.where(pert_gene == d)[0]
        if len(idx) == 0:
            rows.append(dict(dataset=ds, driver=d, in_screen=False)); continue
        o = a.obs.iloc[idx[0]]; z = X.iloc[idx[0]]
        zr = z - np.dot(z - z.mean(), g0) / np.dot(g0, g0) * generic
        pos = z.index.isin(pos_sets[d] - {d}); bg = ~z.index.isin(all_sets[d] | {d})
        med = {k: zr[zr.index.isin(v - {d})].median() for k, v in pos_sets.items()}
        others = [v for k, v in med.items() if k != d]
        rows.append(dict(dataset=ds, driver=d, in_screen=True, n_cells=int(o.num_cells_filtered),
                         kd_fold_expr=round(float(o.fold_expr), 3), cor_with_generic=round(float(np.corrcoef(z, generic)[0, 1]), 3),
                         n_pos_measured=int(pos.sum()), med_z_pos=round(z[pos].median(), 3), med_z_bg=round(z[bg].median(), 3),
                         p_down=mw(z[pos], z[bg]), med_zres_pos=round(zr[pos].median(), 3), med_zres_bg=round(zr[bg].median(), 3),
                         p_down_resid=mw(zr[pos], zr[bg]),
                         own_rank=1 + sum(v < med[d] for v in others), n_regulons=len(others) + 1))
res = pd.DataFrame(rows)
for c in ["p_down", "p_down_resid"]:  # BH within each screen
    m = res[c].notna()
    res.loc[m, c + "_BH"] = res[m].groupby("dataset")[c].transform(lambda p: false_discovery_control(p.values))
res.to_csv(os.path.join(OUT_DIR, "perturbseq_regulon_tests.tsv"), sep="\t", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_rows", 200)
t = res[res.in_screen == True]
print(t[["dataset", "driver", "n_cells", "kd_fold_expr", "cor_with_generic", "n_pos_measured", "med_z_pos", "med_z_bg", "p_down_BH", "med_zres_pos", "med_zres_bg", "p_down_resid_BH", "own_rank", "n_regulons"]].to_string(index=False))
print("\nPer screen: candidates in the screen | p_down_BH<0.05 | p_down_resid_BH<0.05 | own regulon in the first third")
for ds, g in t.groupby("dataset"):
    print(ds, len(g), (g.p_down_BH < 0.05).sum(), (g.p_down_resid_BH < 0.05).sum(), (g.own_rank <= g.n_regulons / 3).sum())
print("not knocked down in each screen:", res[res.in_screen == False].groupby("dataset").driver.apply(lambda s: ", ".join(s)).to_dict())
