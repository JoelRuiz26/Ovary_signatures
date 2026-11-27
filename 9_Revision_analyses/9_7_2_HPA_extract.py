#!/usr/bin/env python
"""Extracts the IHC of normal ovary and ovarian cancer from the Human Protein Atlas (v25) XML entries.

Input:  $OVARY_SIGNATURES/9_Revision_analyses/hpa/xml/<ENSG>.xml (9_7_1_HPA_download_xml.sh),
        9_7_0_HPA_genes.tsv (symbol, Ensembl gene ID; next to this script)
Output ($OVARY_SIGNATURES/9_Revision_analyses/hpa/):
  hpa_v25_ovary_by_antibody.tsv   one row per gene x antibody
  hpa_v25_ovary_patients.tsv      one row per gene x antibody x patient image (normal and tumor)
Usage:  python 9_7_2_HPA_extract.py [xml_dir] [out_dir]
"""
import os, sys, csv, collections
import xml.etree.ElementTree as ET
from pathlib import Path

BASE = Path(os.path.expanduser(os.environ.get("OVARY_SIGNATURES", "~/Ovary_signatures")))  # project data root
HPA_DIR = BASE / "9_Revision_analyses/hpa"
GENES = dict(l.split() for l in open(Path(__file__).parent / "9_7_0_HPA_genes.tsv") if l.strip())  # symbol -> Ensembl
ORD_ST = ["not detected", "low", "medium", "high"]
ORD_IN = ["negative", "weak", "moderate", "strong"]


def txt(el, tag, **attr):
    """Text of the first <tag> child with those attributes (or '')."""
    for c in el.findall(tag):
        if all(c.get(k) == v for k, v in attr.items()):
            return (c.text or "").strip()
    return ""


def count_ordered(vals, order):
    """'not detected 2; low 3' in the given order; unexpected values at the end."""
    c = collections.Counter(v.lower() for v in vals)
    keys = [k for k in order if k in c] + sorted(k for k in c if k not in order)
    return "; ".join(f"{k} {c[k]}" for k in keys)


def count_free(vals):
    c = collections.Counter(vals)
    return "; ".join(f"{k} {n}" for k, n in sorted(c.items(), key=lambda x: (-x[1], x[0])))


def main(xml_dir, out_dir):
    rows, patients = [], []
    for gene, ensg in GENES.items():
        e = ET.parse(Path(xml_dir) / f"{ensg}.xml").getroot().find("entry")
        if e.findtext("name") != gene:
            print("WARNING: different symbol in HPA:", gene, e.findtext("name"), file=sys.stderr)
        ident = e.find("identifier")
        base = {
            "gene": gene, "ensembl": ensg,
            "hpa_entry_version": e.get("version"), "ensembl_version": ident.get("version"),
            "url_gene": e.get("url"),
            "url_cancer_page": f"https://v{e.get('version')}.proteinatlas.org/{ensg}-{gene}/cancer/ovarian+cancer",
        }
        # --- gene level: consensus expression in normal ovary and reliability ---
        g_fol = g_str = g_rel = g_rel_desc = ""
        for te in e.findall("tissueExpression"):
            v = te.find("verification")
            if v is not None:
                g_rel, g_rel_desc = (v.text or "").strip(), (v.get("description") or "").strip()
            for d in te.findall("data"):
                if d.findtext("tissue") == "Ovary":
                    for tc in d.findall("tissueCell"):
                        ct = tc.findtext("cellType")
                        if ct == "Follicle cells":
                            g_fol = txt(tc, "level", type="expression")
                        elif ct == "Ovarian stroma cells":
                            g_str = txt(tc, "level", type="expression")
        base.update(gene_normal_stroma_expression=g_str, gene_normal_follicle_expression=g_fol,
                    gene_ihc_reliability=g_rel, gene_ihc_reliability_note=g_rel_desc)
        # --- prognosis (RNA, TCGA) ---
        prog = {}
        for ce in e.findall("cancerExpression"):
            for d in ce.findall("data"):
                t = d.findtext("tissue") or ""
                if t.startswith("Ovary Serous"):
                    s = d.find("survivalAnalysis")
                    prog[t] = s.attrib if s is not None else {}
        for label, pref in (("Ovary Serous Cystadenocarcinoma (TCGA)", "prog_tcga"),
                               ("Ovary Serous Cystadenocarcinoma (validation)", "prog_validation")):
            a = prog.get(label, {})
            base[pref] = a.get("prognostic", "")
            base[pref + "_type"] = a.get("prognosticType", "")
            base[pref + "_p"] = a.get("pValue", "")
        base["prog_other_ovary"] = "; ".join(k for k in prog if k not in (
            "Ovary Serous Cystadenocarcinoma (TCGA)", "Ovary Serous Cystadenocarcinoma (validation)"))

        abs_ = e.findall("antibody")
        if not abs_:
            rows.append({**base, "antibody": "", "note": "no antibodies in the XML"})
        for ab in abs_:
            f = {**base, "antibody": ab.get("id")}
            has_ihc_tissue = has_ihc_cancer = False
            for te in ab.findall("tissueExpression"):
                assay = te.get("assayType")
                if assay == "tissue":
                    has_ihc_tissue = True
                    v = te.find("verification")
                    f["ab_ihc_validation"] = (v.text or "").strip() if v is not None else ""
                    for d in te.findall("data"):
                        if d.findtext("tissue") != "Ovary":
                            continue
                        for tc in d.findall("tissueCell"):
                            ct = tc.findtext("cellType")
                            pref = {"Follicle cells": "normal_follicle",
                                    "Ovarian stroma cells": "normal_stroma"}.get(ct)
                            if pref is None:
                                f["normal_other_cell_types"] = (f.get("normal_other_cell_types", "") + ct + ";")
                                continue
                            f[pref + "_staining"] = txt(tc, "level", type="staining")
                            f[pref + "_intensity"] = txt(tc, "level", type="intensity")
                            f[pref + "_quantity"] = tc.findtext("quantity") or ""
                            f[pref + "_location"] = tc.findtext("location") or ""
                        pts = d.findall("patient")
                        f["normal_n_patients"] = len(pts)
                        for p in pts:
                            for s in p.findall("sample"):
                                sn = "|".join(x.get("tissueDescription") for x in s.iter("snomed"))
                                for u in s.iter("imageUrl"):
                                    patients.append(dict(
                                        gene=gene, antibody=ab.get("id"), tissue="Ovary (normal)",
                                        patientId=p.findtext("patientId"), sex=p.findtext("sex"),
                                        age=p.findtext("age"), staining="", intensity="", quantity="",
                                        location="", snomed=sn, imageUrl=u.text.strip()))
                elif assay == "cancer":
                    has_ihc_cancer = True
                    for d in te.findall("data"):
                        if d.findtext("tissue") != "Ovarian cancer":
                            continue
                        tc = d.find("tissueCell")
                        summary = {l.text.strip().lower(): int(l.get("count"))
                                   for l in tc.findall("level") if l.get("count")} if tc is not None else {}
                        pts = d.findall("patient")
                        st = [txt(p, "level", type="staining") for p in pts]
                        it = [txt(p, "level", type="intensity") for p in pts]
                        f["tumor_n_patients"] = len(pts)
                        f["tumor_staining"] = count_ordered(st, ORD_ST)
                        f["tumor_intensity"] = count_ordered(it, ORD_IN)
                        f["tumor_quantity"] = count_free([p.findtext("quantity") or "" for p in pts])
                        f["tumor_location"] = count_free([p.findtext("location") or "" for p in pts])
                        c_st = collections.Counter(x.lower() for x in st)
                        c_in = collections.Counter(x.lower() for x in it)
                        f["tumor_n_staining_detected"] = len(pts) - c_st["not detected"]
                        f["tumor_n_staining_med_high"] = c_st["medium"] + c_st["high"]
                        f["tumor_n_intensity_nonneg"] = len(pts) - c_in["negative"]
                        # internal check: summary in the XML vs. sum over patients
                        f["xml_summary_matches"] = str(summary == dict(c_st))
                        hist = []
                        for p in pts:
                            sns = set()
                            urls = []
                            for s in p.findall("sample"):
                                sns.add("|".join(x.get("tissueDescription") for x in s.iter("snomed")))
                                urls += [u.text.strip() for u in s.iter("imageUrl")]
                            h = " // ".join(sorted(sns))
                            hist.append(h)
                            for u in urls:
                                patients.append(dict(
                                    gene=gene, antibody=ab.get("id"), tissue="Ovarian cancer",
                                    patientId=p.findtext("patientId"), sex=p.findtext("sex"),
                                    age=p.findtext("age"), staining=txt(p, "level", type="staining"),
                                    intensity=txt(p, "level", type="intensity"),
                                    quantity=p.findtext("quantity") or "",
                                    location=p.findtext("location") or "", snomed=h, imageUrl=u))
                        f["tumor_histology"] = count_free(hist)
            f["ab_has_ihc_tissue"] = has_ihc_tissue
            f["ab_has_ihc_cancer"] = has_ihc_cancer
            rows.append(f)

    cols = ["gene", "ensembl", "antibody", "ab_has_ihc_tissue", "ab_has_ihc_cancer", "ab_ihc_validation",
            "normal_n_patients",
            "normal_stroma_staining", "normal_stroma_intensity", "normal_stroma_quantity", "normal_stroma_location",
            "normal_follicle_staining", "normal_follicle_intensity", "normal_follicle_quantity",
            "normal_follicle_location", "normal_other_cell_types",
            "tumor_n_patients", "tumor_staining", "tumor_intensity", "tumor_quantity", "tumor_location",
            "tumor_n_staining_detected", "tumor_n_staining_med_high", "tumor_n_intensity_nonneg",
            "tumor_histology", "xml_summary_matches",
            "gene_normal_stroma_expression", "gene_normal_follicle_expression",
            "gene_ihc_reliability", "gene_ihc_reliability_note",
            "prog_tcga", "prog_tcga_type", "prog_tcga_p",
            "prog_validation", "prog_validation_type", "prog_validation_p", "prog_other_ovary",
            "hpa_entry_version", "ensembl_version", "url_gene", "url_cancer_page", "note"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ver = rows[0]["hpa_entry_version"]  # major version declared in the XML
    with open(out / f"hpa_v{ver}_ovary_by_antibody.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval="NA")
        w.writeheader()
        for f in rows:
            w.writerow({k: ("NA" if v == "" else v) for k, v in f.items()})
    pc = ["gene", "antibody", "tissue", "patientId", "sex", "age", "staining", "intensity",
          "quantity", "location", "snomed", "imageUrl"]
    with open(out / f"hpa_v{ver}_ovary_patients.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, pc, delimiter="\t")
        w.writeheader()
        w.writerows(patients)
    print(f"{len(rows)} gene x antibody rows; {len(patients)} patient-image rows")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else HPA_DIR / "xml",
         sys.argv[2] if len(sys.argv) > 2 else HPA_DIR)
