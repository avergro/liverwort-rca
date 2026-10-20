#!/usr/bin/env python3
"""KEGG module completeness and KO breadth per species, from a uniform KOfam annotation.

Asks whether liverworts hold the metabolic functions of angiosperms with fewer
copies, measured on functions (KEGG orthologs, KOs) and on pathways (KEGG
modules) instead of orthogroups, and with an annotation that is not transferred
from Arabidopsis: every proteome is searched with the same KOfam profiles and
the KOfam score thresholds.

Steps
  prepare  parse the module definitions of the snapshot, keep pathway modules,
           write the list of their KOs and extract those profiles (hmmfetch)
  search   hmmsearch of the module-KO profiles against the 11 proteomes
  analyse  assign KOs (KOfam thresholds), evaluate every module per species,
           and summarise breadth, copies per KO and complete modules

Inputs
  data/kofam/KEGG/            anvi'o KEGG snapshot (KEGG_build_2025-02-04): HMMs, ko_list.txt, modules/
  data/proteomes/filtered/    the corrected proteomes used for OrthoFinder (one protein per gene)
Outputs
  data/results_kegg_modules/

Module logic: the definition is split into steps (space-separated at depth 0);
inside a step "," is OR, a space inside parentheses is AND, "+" joins required
complex components and "-" marks an optional component. A step without a KO
("--") is ignored. A module is complete when every step is satisfied, and
near-complete when at most one step is missing (the KEGG criterion for
incomplete genomes).

Usage
  python3 33_kegg_module_completeness.py prepare
  python3 33_kegg_module_completeness.py search --jobs 4 --cpu 3
  python3 33_kegg_module_completeness.py analyse
"""
import argparse
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
KEGG = ROOT / "data" / "kofam" / "KEGG"
PROT = ROOT / "data" / "proteomes" / "filtered"
OUT = ROOT / "data" / "results_kegg_modules"
HMMBIN = Path(os.environ.get("HMMER_BIN", ""))
SPECIES = ["Afiliculoides", "Athaliana", "Creinhardtii", "Lcruciata", "Mpolymorpha", "Pabies",
           "Pamericana", "Ppatens", "Scucullata", "Slycopersicum", "Smoellendorffii"]
LIVER = ["Mpolymorpha", "Lcruciata"]
ANGIO = ["Athaliana", "Slycopersicum", "Pamericana"]
LAND = [s for s in SPECIES if s != "Creinhardtii"]


# ── module definitions ─────────────────────────────────────────────────────

def read_modules():
    """Return {module: dict(name, cls, definition)} for pathway modules."""
    mods = {}
    for f in sorted((KEGG / "modules").glob("M*")):
        txt = f.read_text().splitlines()
        rec, field = {"definition": []}, None
        for line in txt:
            key, val = line[:12].strip(), line[12:].strip()
            if key:
                field = key
            if field == "NAME" and key:
                rec["name"] = val
            elif field == "CLASS" and key:
                rec["cls"] = val
            elif field == "DEFINITION":
                rec["definition"].append(val)
        if rec.get("cls", "").startswith("Pathway modules"):
            rec["definition"] = " ".join(rec["definition"])
            mods[f.name] = rec
    return mods


def tokenize(s):
    return re.findall(r"K\d{5}|M\d{5}|--|[(),+\- ]", s)


def split_steps(defn):
    """Split a definition at depth-0 spaces."""
    steps, depth, cur = [], 0, ""
    for ch in defn:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == " " and depth == 0:
            if cur:
                steps.append(cur)
            cur = ""
        else:
            cur += ch
    if cur:
        steps.append(cur)
    return steps


class Expr:
    """Recursive-descent evaluator of one step against a set of KOs (and complete modules)."""

    def __init__(self, text, present):
        self.t = [x for x in tokenize(text)]
        self.i, self.present = 0, present

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None

    def take(self):
        self.i += 1
        return self.t[self.i - 1]

    def expr(self):                     # OR of terms
        v = self.term()
        while self.peek() == ",":
            self.take()
            v = self.term() or v
        return v

    def term(self):                     # AND of factors separated by spaces
        v = self.factor()
        while self.peek() == " ":
            self.take()
            if self.peek() in (")", ",", None):
                break
            v = self.factor() and v
        return v

    def factor(self):                   # units joined by + (required) or - (optional)
        v = self.unit()
        while self.peek() in ("+", "-"):
            op = self.take()
            u = self.unit()
            if op == "+":
                v = v and u
        return v

    def unit(self):
        tok = self.take()
        if tok == "(":
            v = self.expr()
            assert self.take() == ")"
            return v
        if tok == "--":
            return True
        if tok == "-":                  # leading optional component
            self.unit()
            return True
        return tok in self.present


def evaluate(defn, present):
    """Return (n_steps, n_satisfied) ignoring KO-less steps."""
    n = ok = 0
    for st in split_steps(defn):
        if st == "--":
            continue
        n += 1
        ok += bool(Expr(st, present).expr())
    return n, ok


# ── steps ──────────────────────────────────────────────────────────────────

def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    mods = read_modules()
    kos = sorted({k for m in mods.values() for k in re.findall(r"K\d{5}", m["definition"])})
    kl = pd.read_csv(KEGG / "ko_list.txt", sep="\t")
    hmm = sorted((KEGG / "HMMs").glob("*.hmm"))[0]
    in_hmm = set(subprocess.run(["grep", "^NAME", hmm], capture_output=True, text=True).stdout.split()[1::2])
    have = set(kl.knum[kl.threshold != "-"]) & in_hmm
    kos_ok = [k for k in kos if k in have]
    missing = [k for k in kos if k not in have]
    (OUT / "module_kos.txt").write_text("\n".join(kos_ok) + "\n")
    (OUT / "module_kos_without_profile.txt").write_text("\n".join(missing) + "\n")
    pd.DataFrame([(m, r["name"], r["cls"], r["definition"]) for m, r in mods.items()],
                 columns=["module", "name", "class", "definition"]).to_csv(OUT / "modules.tsv", sep="\t", index=False)
    if not Path(str(hmm) + ".ssi").exists():
        subprocess.run([HMMBIN / "hmmfetch", "--index", hmm], check=True)
    with open(OUT / "module_kos.hmm", "w") as fh:
        subprocess.run([HMMBIN / "hmmfetch", "-f", hmm, OUT / "module_kos.txt"], stdout=fh, check=True)
    print(f"{len(mods)} pathway modules, {len(kos)} KOs, {len(kos_ok)} with a KOfam profile -> {OUT / 'module_kos.hmm'}")


def search(jobs, cpu):
    (OUT / "hmmsearch").mkdir(parents=True, exist_ok=True)

    def run(sp):
        out = OUT / "hmmsearch" / f"{sp}.tbl"
        if out.exists() and out.read_text().rstrip().endswith("[ok]"):
            return f"{sp} (cached)"
        subprocess.run([HMMBIN / "hmmsearch", "--cpu", str(cpu), "--noali", "-E", "1e-3",
                        "--tblout", out, "-o", "/dev/null", OUT / "module_kos.hmm",
                        PROT / f"{sp}.primary.fa"], check=True)
        return sp

    with ThreadPoolExecutor(jobs) as ex:
        for r in ex.map(run, SPECIES):
            print("done", r, flush=True)


def assign():
    kl = pd.read_csv(KEGG / "ko_list.txt", sep="\t", usecols=["knum", "threshold", "score_type"])
    kl = kl[kl.threshold != "-"].astype({"threshold": float}).set_index("knum")
    rows = []
    for sp in SPECIES:
        for line in open(OUT / "hmmsearch" / f"{sp}.tbl"):
            if line.startswith("#"):
                continue
            f = line.split()
            gene, ko, full, dom = f[0], f[2], float(f[5]), float(f[8])
            if ko not in kl.index:
                continue
            thr, st = kl.at[ko, "threshold"], kl.at[ko, "score_type"]
            if (full if st == "full" else dom) >= thr:
                rows.append((sp, gene, ko))
    a = pd.DataFrame(rows, columns=["species", "gene", "ko"]).drop_duplicates()
    a.to_csv(OUT / "ko_assignments.tsv", sep="\t", index=False)
    return a


def analyse():
    a = assign()
    mods = pd.read_csv(OUT / "modules.tsv", sep="\t")
    pres = {sp: set(a.ko[a.species == sp]) for sp in SPECIES}
    # modules defined from other modules (metabolic capacity): add complete modules and re-evaluate
    for _ in range(3):
        for sp in SPECIES:
            done = {m.module for _, m in mods.iterrows()
                    if (lambda r: r[0] == r[1])(evaluate(m.definition, pres[sp]))}
            pres[sp] |= done
    ev = []
    for _, m in mods.iterrows():
        for sp in SPECIES:
            n, ok = evaluate(m.definition, pres[sp])
            ev.append((m.module, m["name"], m["class"].split("; ")[1], sp, n, ok))
    ev = pd.DataFrame(ev, columns=["module", "name", "category", "species", "steps", "satisfied"])
    ev["complete"] = ev.satisfied == ev.steps
    ev["near_complete"] = ev.satisfied >= ev.steps - 1
    ev["fraction"] = ev.satisfied / ev.steps
    ev.to_csv(OUT / "module_completeness.tsv", sep="\t", index=False)

    # plant-relevant modules: complete in at least one land plant
    wide = ev.pivot(index="module", columns="species", values="complete")
    relevant = wide.index[wide[LAND].any(axis=1)]
    kos_rel = {k for m in mods[mods.module.isin(relevant)].definition for k in re.findall(r"K\d{5}", m)}
    rows = []
    for sp in SPECIES:
        g = a[(a.species == sp) & a.ko.isin(kos_rel)]
        e = ev[(ev.species == sp) & ev.module.isin(relevant)]
        rows.append(dict(species=sp, kos=g.ko.nunique(), genes=g.gene.nunique(),
                         genes_per_ko=g.groupby("ko").gene.nunique().mean(),
                         complete=int(e.complete.sum()), near_complete=int(e.near_complete.sum()),
                         mean_fraction=e.fraction.mean()))
    summ = pd.DataFrame(rows).set_index("species")
    summ.to_csv(OUT / "species_summary.tsv", sep="\t", float_format="%.4g")
    ratio = summ.loc[LIVER].mean() / summ.loc[ANGIO].mean()
    print(f"{len(relevant)} modules complete in at least one land plant; {len(kos_rel)} KOs in them\n")
    print(summ.round(3).to_string())
    print("\nliverwort / angiosperm mean:\n" + ratio.round(3).to_string())

    # modules that separate the groups
    sep = wide.loc[relevant]
    only_angio = sep.index[sep[ANGIO].all(axis=1) & ~sep[LIVER].any(axis=1)]
    only_liver = sep.index[sep[LIVER].all(axis=1) & ~sep[ANGIO].any(axis=1)]
    names = mods.set_index("module")
    diff = pd.DataFrame([(m, "complete in all angiosperms, in no liverwort", names.at[m, "name"]) for m in only_angio]
                        + [(m, "complete in both liverworts, in no angiosperm", names.at[m, "name"]) for m in only_liver],
                        columns=["module", "pattern", "name"])
    nfrac = ev.pivot(index="module", columns="species", values="fraction").loc[diff.module, LIVER + ANGIO].round(2)
    diff = diff.join(nfrac, on="module")
    diff.to_csv(OUT / "modules_differing.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 70):
        print("\n" + diff.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "search", "analyse"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpu", type=int, default=3)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"prepare": prepare, "search": lambda: search(a.jobs, a.cpu), "analyse": analyse}[a.step]()
