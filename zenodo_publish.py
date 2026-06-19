#!/usr/bin/env python3
"""
zenodo_publish.py
=================
Archive this repository to Zenodo and obtain a citable DOI.

Requirements
------------
- Python 3.8+ (standard library only, no pip installs needed)
- git installed and working in this directory
- A Zenodo personal access token:
    1. Log in at https://zenodo.org  (use your GitHub or ORCID account)
    2. Go to https://zenodo.org/account/settings/applications/
    3. Click "New token", give it any name, select scopes:
         deposit:write   deposit:actions
    4. Copy the token — it is shown only once.
- ZENODO_TOKEN environment variable set to that token.

Usage
-----
# Test first with the sandbox (no real DOI, no cost):
    ZENODO_TOKEN=<sandbox-token> python3 zenodo_publish.py --sandbox

# Production run (issues a real, permanent DOI):
    ZENODO_TOKEN=<production-token> python3 zenodo_publish.py

Notes
-----
- The sandbox and production sites use separate accounts/tokens.
  Get a sandbox token at https://sandbox.zenodo.org/account/settings/applications/
- Fill in ALL authors in the METADATA block below before the production run.
- The script creates a git tag, pushes it, archives the repo, uploads to
  Zenodo, sets metadata, and publishes — printing the DOI at the end.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  CONFIGURATION — edit this section before running
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

GIT_TAG     = "v1.0.0"
GITHUB_REPO = "JoelRuiz26/Ovary_signatures"
GITHUB_URL  = f"https://github.com/{GITHUB_REPO}"

METADATA: dict = {
    "title": (
        "Integrative Transcriptomics Addresses Control Ambiguity to Identify "
        "a Robust Ovarian Cancer Signature and Functionally Essential Driver "
        "Regulators: analysis code"
    ),
    "upload_type": "software",
    "description": (
        "<p>R analysis pipeline supporting the manuscript submitted to "
        "<em>NAR Cancer &amp; Computational Biology</em> (Oxford Academic). "
        "The workflow integrates multi-cohort transcriptomics (GEO microarray "
        "studies and TCGA RNA-seq), network-based Master Regulator Analysis "
        "(MRA) applied to the full space of transcriptional modulators "
        "(sequence-specific TFs, cofactors, chromatin remodelers, and "
        "epigenetic regulators), and CRISPR-Cas9 dependency data from the "
        "Cancer Dependency Map (DepMap) to derive a robust ovarian cancer "
        "gene signature and a high-confidence set of upstream driver "
        "regulators. The framework identified 17 central transcriptional "
        "modulators that are both master regulators of essential gene programs "
        "and functionally required in ovarian cancer cell lines.</p>"
    ),
    "license": "gpl-3.0",
    "access_right": "open",
    "language": "eng",
    "keywords": [
        "ovarian cancer",
        "transcriptomics",
        "master regulator analysis",
        "DepMap",
        "CRISPR dependency",
        "gene signature",
        "meta-analysis",
        "R bioinformatics",
    ],
    # ── AUTHORS ───────────────────────────────────────────────────────────────
    # Add every co-author before the production run.
    # Format: {"name": "Apellido, Nombre", "affiliation": "...", "orcid": "0000-..."}
    # ORCID is optional but strongly recommended by NAR.
    "creators": [
        {"name": "Ruiz-Hernández, Joel",
         "affiliation": "Programa de Maestría y Doctorado en Ciencias Bioquímicas, UNAM; Instituto Nacional de Medicina Genómica, Mexico City, Mexico",
         "orcid": "0009-0007-8914-6100"},
        {"name": "Pérez-Calixto, Daniel",
         "affiliation": "Departamento de Física, Facultad de Ciencias, UNAM; Instituto Nacional de Medicina Genómica, Mexico City, Mexico",
         "orcid": "0009-0004-3988-1260"},
        {"name": "Pastelín-Morales, José Manuel",
         "affiliation": "Maestría en Biotecnología, Universidad del Papaloapan, Tuxtepec, Oaxaca, Mexico",
         "orcid": "0009-0004-9461-9898"},
        {"name": "Hernández-Lemus, Enrique",
         "affiliation": "Instituto Nacional de Medicina Genómica, Mexico City, Mexico",
         "orcid": "0000-0002-1872-1397"},
        {"name": "Vazquez-Victorio, Genaro",
         "affiliation": "Departamento de Física, Facultad de Ciencias, UNAM; Laboratorio Nacional de Soluciones Biomiméticas, UNAM, Mexico City, Mexico"},
        {"name": "Martínez-Ramírez, Angélica S.",
         "affiliation": "Instituto de Biotecnología, Universidad del Papaloapan, Tuxtepec, Oaxaca, Mexico",
         "orcid": "0000-0002-2474-4739"},
        {"name": "Tovar, Hugo",
         "affiliation": "Departamento de Física, Facultad de Ciencias, UNAM; Instituto Nacional de Medicina Genómica, Mexico City, Mexico",
         "orcid": "0000-0002-8360-6133"},
    ],
    "related_identifiers": [
        {
            "relation":   "isSupplementTo",
            "identifier": GITHUB_URL,
            "scheme":     "url",
        }
    ],
}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  Internal helpers — nothing to edit below this line
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SANDBOX  = "--sandbox" in sys.argv
BASE_URL = "https://sandbox.zenodo.org" if SANDBOX else "https://zenodo.org"

TOKEN = os.environ.get("ZENODO_TOKEN", "")
if not TOKEN:
    base = BASE_URL
    print("Error: ZENODO_TOKEN is not set.\n")
    print(f"  1. Log in at {base}")
    print(f"  2. Go to   {base}/account/settings/applications/")
    print("  3. Click 'New token' — scopes needed: deposit:write  deposit:actions")
    print("  4. Run:  ZENODO_TOKEN=<token> python3 zenodo_publish.py")
    sys.exit(1)


def _http(method: str, url: str, *, json_body=None, binary_body: bytes | None = None) -> dict:
    headers: dict[str, str] = {"Authorization": f"Bearer {TOKEN}"}
    body: bytes | None = None
    if binary_body is not None:
        headers["Content-Type"] = "application/octet-stream"
        body = binary_body
    elif json_body is not None:
        headers["Content-Type"] = "application/json"
        body = json.dumps(json_body).encode()
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")
        print(f"\n✗ HTTP {exc.code} on {method} {url}")
        try:
            parsed = json.loads(detail)
            print(json.dumps(parsed, indent=2))
        except json.JSONDecodeError:
            print(detail)
        sys.exit(1)


def zenodo(method: str, path: str, **kwargs) -> dict:
    return _http(method, f"{BASE_URL}/api{path}", **kwargs)


def git(*args: str, capture: bool = False) -> str:
    if capture:
        return subprocess.check_output(["git", *args], stderr=subprocess.DEVNULL).decode().strip()
    subprocess.check_call(["git", *args])
    return ""


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  Main flow
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if SANDBOX:
    print("⚠  SANDBOX mode — no permanent DOI will be issued\n")

# 1. Tag and push ──────────────────────────────────────────────────────────────
print(f"[1/5] Tagging HEAD as {GIT_TAG}...")
existing_tags = git("tag", "--list", GIT_TAG, capture=True)
if existing_tags:
    print(f"      Tag {GIT_TAG} already exists — skipping creation.")
else:
    git("tag", GIT_TAG)
git("push", "origin", GIT_TAG)
print(f"      Tag {GIT_TAG} pushed to {GITHUB_URL}.")

# 2. Archive ───────────────────────────────────────────────────────────────────
archive_name = f"ovary_signatures_{GIT_TAG}.tar.gz"
archive_path = Path("/tmp") / archive_name
print(f"\n[2/5] Archiving repository ({GIT_TAG}) → {archive_path}...")
git("archive", "--format=tar.gz", f"--output={archive_path}", GIT_TAG)
size_kb = archive_path.stat().st_size // 1024
print(f"      Archive ready: {size_kb} KB")

# 3. Create deposition ─────────────────────────────────────────────────────────
print("\n[3/5] Creating Zenodo deposition...")
dep     = zenodo("POST", "/deposit/depositions", json_body={})
dep_id  = dep["id"]
bucket  = dep["links"]["bucket"]
print(f"      Deposition ID: {dep_id}")

# 4. Upload archive ────────────────────────────────────────────────────────────
print(f"\n[4/5] Uploading {archive_name} to Zenodo...")
raw = archive_path.read_bytes()
uploaded = _http("PUT", f"{bucket}/{archive_name}", binary_body=raw)
print(f"      Uploaded: {uploaded['key']}  ({uploaded['size'] // 1024} KB)")
archive_path.unlink()

# 5. Set metadata ──────────────────────────────────────────────────────────────
print("\n[5/5] Setting metadata and publishing...")
zenodo("PUT", f"/deposit/depositions/{dep_id}", json_body={"metadata": METADATA})

published   = zenodo("POST", f"/deposit/depositions/{dep_id}/actions/publish")
doi         = published["doi"]
concept_doi = published.get("conceptdoi", doi)
record_url  = published["links"]["html"]

# ── Report ────────────────────────────────────────────────────────────────────
year  = date.today().year
title = METADATA["title"]
author_block = ", ".join(
    c["name"].split(",")[0] for c in METADATA["creators"]
)

print()
print("━" * 68)
print("  SUCCESS" + ("  (sandbox)" if SANDBOX else ""))
print("━" * 68)
print(f"  Version DOI : https://doi.org/{doi}")
print(f"  Concept DOI : https://doi.org/{concept_doi}")
print(f"               └─ cite this one in the paper (tracks all versions)")
print(f"  Zenodo page : {record_url}")
print()
print("  NAR reference entry:")
print(f"  {author_block} et al. ({year}) {title}.")
print(f"  Zenodo. https://doi.org/{concept_doi}.")
print()
print("  Code Availability sentence:")
print(f"  All custom R scripts are publicly available on GitHub")
print(f"  ({GITHUB_URL}) and archived at Zenodo")
print(f"  (https://doi.org/{concept_doi}).")
print("━" * 68)
print()
print("  Next steps:")
print("  1. Open the Zenodo page above and verify authors / metadata.")
print("  2. Add the concept DOI to the Code Availability section.")
print("  3. Add the reference to your NAR reference list.")
