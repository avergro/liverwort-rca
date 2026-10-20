#!/usr/bin/env python3
"""
BUSCO completeness of every proteome, in one table read by the later scripts.

Reads the short_summary*.json written by 03_busco_qc.sh under data/busco/ and
writes data/busco/busco_summary.tsv (species code, species name, lineage
dataset and the C, S, D, F and M percentages). Chlamydomonas was assessed with
both lineage datasets; the chlorophyta run is kept, which is the applicable one.

Usage:
    python3 03b_busco_summary.py
"""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
BUSCO = ROOT / "data" / "busco"
NAMES = {"Mpolymorpha": "Marchantia polymorpha", "Lcruciata": "Lunularia cruciata",
         "Smoellendorffii": "Selaginella moellendorffii", "Ppatens": "Physcomitrium patens",
         "Afiliculoides": "Azolla filiculoides", "Scucullata": "Salvinia cucullata",
         "Creinhardtii": "Chlamydomonas reinhardtii", "Athaliana": "Arabidopsis thaliana",
         "Slycopersicum": "Solanum lycopersicum", "Pamericana": "Persea americana",
         "Pabies": "Picea abies"}


def main():
    rows = []
    for d in sorted(p for p in BUSCO.iterdir() if p.is_dir()):
        js = list(d.glob("short_summary*.json"))
        if not js:
            continue
        data = json.loads(js[0].read_text())
        res = data.get("results", {})
        code = d.name.replace("_chlorophyta", "")
        rows.append({"species": code, "Species": NAMES.get(code, code), "run": d.name,
                     "BUSCO dataset": data.get("lineage_dataset", {}).get("name", ""),
                     "C (%)": res.get("Complete percentage"), "S (%)": res.get("Single copy percentage"),
                     "D (%)": res.get("Multi copy percentage"), "F (%)": res.get("Fragmented percentage"),
                     "M (%)": res.get("Missing percentage")})
    df = pd.DataFrame(rows)
    df = df.sort_values("run", key=lambda c: c.str.contains("chlorophyta"))
    df = df.drop_duplicates("species", keep="last").drop(columns="run").sort_values("species")
    df.to_csv(BUSCO / "busco_summary.tsv", sep="\t", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
