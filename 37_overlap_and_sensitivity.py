#!/usr/bin/env python3
"""
Orthogroup overlap with the angiosperms, baselines and sensitivity checks.

Coverage parity says that liverworts occupy as many orthogroups as the
angiosperms, not that they occupy the same ones. This script reports which
orthogroups each species shares with the angiosperms, the baseline against
which the Filter B sharing with the moss and the lycophyte should be read, three
sensitivity checks of the coverage and redundancy ratios (hierarchical
orthogroups at the land plant node, fragmentary genes removed, lineage-specific
orthogroups removed), the gene-model counts of the Lunularia input, and the
correlation that fixes the sign of the presence ranking.

Outputs (data/results_controls/):
    M_overlap_with_angiosperms.tsv
    M_filterB_sharing_baseline.tsv
    M_sensitivity_ratios.tsv
    M_lunularia_input.tsv
    M_sign_convention.tsv
    M_lineage_losses.tsv, M_deficit_decomposition.tsv
    M_regulatory_stratified.tsv (copy ratio of regulatory and metabolic orthogroups by angiosperm copies)
    M_missing_GO_enrichment.tsv (needs the GOA and go-basic files in CODE/data/go)
"""

import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OF = ROOT / "data/orthofinder_v3/Results_Oct05"
# hierarchical orthogroups from gene trees reconciled with the accepted species tree (38_fixed_species_tree.py)
TREES = ROOT / "data/orthofinder_v3/Results_Oct05_fixedtree"
OUT = ROOT / "data/results_controls"
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]
OTHER_NON_SEED = ["Ppatens", "Smoellendorffii", "Afiliculoides", "Scucullata"]
LINEAGES = {"liverworts": LIVERWORTS, "angiosperms": ANGIOSPERMS,
            "ferns": ["Afiliculoides", "Scucullata"]}
CONTEXTS = ["primary_metabolism", "C_metabolism", "N_metabolism"]


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ratios(counts):
    """Liverwort / angiosperm ratios of orthogroups occupied and genes per occupied orthogroup."""
    occ = (counts > 0).sum()
    gpo = counts.sum() / occ
    return (occ[LIVERWORTS].mean() / occ[ANGIOSPERMS].mean(),
            gpo[LIVERWORTS].mean() / gpo[ANGIOSPERMS].mean())


def overlap(counts):
    P = counts > 0
    all_ang = P[ANGIOSPERMS].all(axis=1)
    any_ang = P[ANGIOSPERMS].any(axis=1)
    rows = []
    for sp in counts.columns:
        # orthogroups in all three angiosperms and in at least one other non-seed land plant
        others = [s for s in OTHER_NON_SEED + LIVERWORTS if s != sp]
        conserved = all_ang & P[others].any(axis=1)
        # reference that treats every lineage alike: orthogroups present in every land plant
        # outside the lineage of the species (so losses shared within a lineage are seen)
        in_all_others = P[outside_lineage(sp, counts.columns)].all(axis=1)
        rows.append({"species": sp, "occupied": int(P[sp].sum()),
                     "in_all_three_angiosperms": int(all_ang.sum()),
                     "missing_of_all_three_angiosperms": int((all_ang & ~P[sp]).sum()),
                     "conserved_with_non_seed": int(conserved.sum()),
                     "missing_of_conserved": int((conserved & ~P[sp]).sum()),
                     "absent_from_all_angiosperms": int((P[sp] & ~any_ang).sum()),
                     "in_all_land_plants_outside_lineage": int(in_all_others.sum()),
                     "missing_of_those_outside_lineage": int((in_all_others & ~P[sp]).sum()),
                     "jaccard_with_Athaliana": (P[sp] & P.Athaliana).sum() / (P[sp] | P.Athaliana).sum()})
    d = pd.DataFrame(rows)
    d["fraction_missing_of_all_three"] = d.missing_of_all_three_angiosperms / d.in_all_three_angiosperms
    d.loc[d.species.isin(ANGIOSPERMS), ["missing_of_all_three_angiosperms", "fraction_missing_of_all_three",
                                        "absent_from_all_angiosperms", "conserved_with_non_seed",
                                        "missing_of_conserved"]] = np.nan
    d["fraction_missing_of_conserved"] = d.missing_of_conserved / d.conserved_with_non_seed
    d["fraction_missing_outside_lineage"] = d.missing_of_those_outside_lineage / d.in_all_land_plants_outside_lineage
    return d


def busco_expected(ov, losses):
    """Loss expected from proteome completeness: linear regression of the fraction of the reference
    missing on BUSCO completeness across the land plants outside the liverworts, for orthogroups and KOs."""
    b = pd.read_csv(ROOT / "data/busco/busco_summary.tsv", sep="\t")      # 03b_busco_summary.py
    busco = b.set_index("species")["C (%)"]
    ko = losses[losses.unit == "KEGG orthologs"]
    ko_frac = {}
    for _, r in ko.iterrows():
        for c in r.index:
            if c.startswith("missing_") and c != "missing_in_whole_lineage" and pd.notna(r[c]):
                ko_frac.setdefault(c[8:], r[c] / r.reference)     # first (own-lineage) row of each species
    rows = []
    og = ov.set_index("species").fraction_missing_outside_lineage
    land = [x for x in og.index if x != "Creinhardtii"]
    fit_on = [x for x in land if x not in LIVERWORTS]
    for unit, y in (("orthogroups", og), ("KEGG orthologs", pd.Series(ko_frac))):
        slope, icpt = np.polyfit(busco[fit_on], y[fit_on], 1)
        resid = y[fit_on] - (icpt + slope * busco[fit_on])
        sd = resid.std(ddof=2)
        for x in land:
            exp = icpt + slope * busco[x]
            rows.append({"unit": unit, "species": x, "busco": busco[x], "observed": y[x], "expected": exp,
                         "residual_sd": (y[x] - exp) / sd, "slope_per_busco_point": slope})
    return pd.DataFrame(rows)


def kegg_best_hit():
    """Module completeness when each gene is assigned only to its best-scoring KO (score / threshold)."""
    spec = importlib.util.spec_from_file_location("m33", ROOT / "33_kegg_module_completeness.py")
    m33 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m33)
    if not (m33.OUT / "hmmsearch").exists():
        return None
    kl = pd.read_csv(m33.KEGG / "ko_list.txt", sep="\t", usecols=["knum", "threshold", "score_type"])
    kl = kl[kl.threshold != "-"].astype({"threshold": float}).set_index("knum")
    rows = []
    for sp in m33.SPECIES:
        for line in open(m33.OUT / "hmmsearch" / f"{sp}.tbl"):
            if line.startswith("#"):
                continue
            f = line.split()
            gene, ko, full, dom = f[0], f[2], float(f[5]), float(f[8])
            if ko in kl.index:
                thr, st = kl.at[ko, "threshold"], kl.at[ko, "score_type"]
                sc = full if st == "full" else dom
                if sc >= thr:
                    rows.append((sp, gene, ko, sc / thr))
    a = pd.DataFrame(rows, columns=["species", "gene", "ko", "r"]).sort_values("r", ascending=False)
    mods = pd.read_csv(m33.OUT / "modules.tsv", sep="\t")
    out = []
    for label, A in (("all KOs above threshold", a.drop_duplicates(["species", "gene", "ko"])),
                     ("best-scoring KO per gene", a.drop_duplicates(["species", "gene"]))):
        pres = {sp: set(A.ko[A.species == sp]) for sp in m33.SPECIES}
        for _ in range(3):
            for sp in m33.SPECIES:
                pres[sp] |= {m.module for _, m in mods.iterrows()
                             if (lambda r: r[0] == r[1])(m33.evaluate(m.definition, pres[sp]))}
        W = pd.DataFrame({sp: {m.module: (lambda r: r[0] == r[1])(m33.evaluate(m.definition, pres[sp]))
                               for _, m in mods.iterrows()} for sp in m33.SPECIES})
        W = W[W[m33.LAND].any(axis=1)]
        a3 = W[m33.ANGIO].all(axis=1)
        gaps = sorted(W.index[a3 & ~W.Mpolymorpha])
        out.append({"assignment": label, "modules_complete_in_a_land_plant": len(W),
                    "complete_in_all_three_angiosperms": int(a3.sum()),
                    "of_those_complete_in_Mpolymorpha": int((a3 & W.Mpolymorpha).sum()),
                    "of_those_complete_in_Lcruciata": int((a3 & W.Lcruciata).sum()),
                    "Mpolymorpha_gaps": ", ".join(gaps)})
    return pd.DataFrame(out)


def outside_lineage(sp, columns):
    own = next((g for g in LINEAGES.values() if sp in g), [sp])
    return [c for c in columns if c not in own and c != "Creinhardtii"]


def lineage_losses(counts):
    """Orthogroups, KEGG orthologs and complete modules present in every land plant outside a
    lineage and absent from all members of that lineage (and from each member)."""
    P = counts > 0
    ko = pd.read_csv(ROOT / "data/results_kegg_modules/ko_assignments.tsv", sep="\t")
    mods = pd.read_csv(ROOT / "data/results_kegg_modules/modules.tsv", sep="\t")
    comp = pd.read_csv(ROOT / "data/results_kegg_modules/module_completeness.tsv", sep="\t")
    land = [c for c in counts.columns if c != "Creinhardtii"]
    wide = comp.pivot(index="module", columns="species", values="complete").astype(bool)
    relevant = set(wide.index[wide[land].any(axis=1)])
    kos_rel = {k for m in mods[mods.module.isin(relevant)].definition for k in re.findall(r"K\d{5}", m)}
    K = pd.crosstab(ko.ko, ko.species) > 0
    K = K.reindex(sorted(kos_rel), fill_value=False)
    mod = pd.read_csv(ROOT / "data/results_kegg_modules/module_completeness.tsv", sep="\t")
    C = mod.pivot(index="module", columns="species", values="complete").astype(bool)
    rows = []
    groups = dict(LINEAGES, moss=["Ppatens"], lycophyte=["Smoellendorffii"], gymnosperm=["Pabies"],
                  bryophytes=LIVERWORTS + ["Ppatens"], seed_plants=ANGIOSPERMS + ["Pabies"])
    for lin, members in groups.items():
        out = [c for c in counts.columns if c not in members and c != "Creinhardtii"]
        for unit, M in (("orthogroups", P), ("KEGG orthologs", K), ("complete modules", C)):
            ref = M[out].all(axis=1)
            row = {"lineage": lin, "unit": unit, "reference": int(ref.sum()),
                   "missing_in_whole_lineage": int((ref & ~M[members].any(axis=1)).sum())}
            for m in members:
                row[f"missing_{m}"] = int((ref & ~M[m]).sum())
            rows.append(row)
    return pd.DataFrame(rows)


def deficit_decomposition(counts):
    """Difference in genes between each angiosperm and each liverwort over the orthogroups
    present in at least two species, split into four additive parts."""
    P = counts > 0
    rows = []
    for L in LIVERWORTS:
        nonseed = [s for s in OTHER_NON_SEED + LIVERWORTS if s != L]
        for a in ANGIOSPERMS:
            total = counts[a].sum() - counts[L].sum()
            shared = P[a] & P[L]
            only_a = P[a] & ~P[L]
            lost = only_a & P[nonseed].any(axis=1)
            rows.append({"liverwort": L, "angiosperm": a, "gene_difference": int(total),
                         "copies_in_shared_orthogroups": (counts.loc[shared, a] - counts.loc[shared, L]).sum() / total,
                         "orthogroups_lacking_but_in_other_non_seed": counts.loc[lost, a].sum() / total,
                         "orthogroups_restricted_to_seed_plants": counts.loc[only_a & ~lost, a].sum() / total,
                         "orthogroups_of_the_liverwort_absent_from_angiosperm": -counts.loc[P[L] & ~P[a], L].sum() / total})
    return pd.DataFrame(rows)


def filterB_baseline():
    """Share of all orthogroups of the context enriched in the moss, the lycophyte or either,
    against the share among the Filter B orthogroups."""
    rows = []
    for ctx in CONTEXTS:
        M = pd.read_csv(ROOT / f"data/results_{ctx}_v3/matrix_rca.tsv", sep="\t", index_col=0)
        b = pd.read_csv(ROOT / "gene_lists" / ctx / "02_OGs_both_liverworts.tsv", sep="\t", index_col=0)
        no_ang = M[(M[ANGIOSPERMS] == 0).all(axis=1)]
        for name, sel in (("all orthogroups of the context", M),
                          ("no angiosperm enriched", no_ang), ("Filter B", b)):
            pp, sm = sel.Ppatens == 1, sel.Smoellendorffii == 1
            rows.append({"context": ctx, "set": name, "n": len(sel),
                         "Ppatens": pp.mean(), "Smoellendorffii": sm.mean(), "either": (pp | sm).mean()})
    return pd.DataFrame(rows)


def hog_counts(node):
    h = pd.read_csv(TREES / "Phylogenetic_Hierarchical_Orthogroups" / f"{node}.tsv", sep="\t",
                    index_col=0, dtype=str)
    sp = [c for c in h.columns if c not in ("OG", "Gene Tree Parent Clade")]
    n = h[sp].apply(lambda col: col.fillna("").map(lambda x: 0 if x == "" else len(x.split(", "))))
    return n[(n > 0).sum(axis=1) >= 2]


def fasta_lengths(fa):
    L, k = {}, None
    for line in open(fa):
        if line.startswith(">"):
            k = line[1:].split()[0].replace(":", "_")   # OrthoFinder writes ":" as "_"
            L[k] = 0
        else:
            L[k] += len(line.strip().rstrip("*"))
    return L


def without_fragments(bm, threshold=0.5):
    """Counts after removing genes shorter than `threshold` times the median length of their orthogroup."""
    og = pd.read_csv(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"), sep="\t", index_col=0)
    L = {}
    for s in og.columns:
        L.update(fasta_lengths(ROOT / "data/orthofinder_input" / f"{s}.fa"))
    rows = [(o, s, L[g]) for o, r in og.iterrows() for s in og.columns if isinstance(r[s], str)
            for g in r[s].split(", ")]
    df = pd.DataFrame(rows, columns=["og", "sp", "len"])
    df["rel"] = df["len"] / df.og.map(df.groupby("og")["len"].median())
    kept = df[df.rel >= threshold].groupby(["og", "sp"]).size().unstack(fill_value=0)
    return kept[(kept > 0).sum(axis=1) >= 2], df


def without_lineage_specific(counts):
    """Counts restricted, for each species, to orthogroups shared with a species outside its lineage."""
    P = counts > 0
    out = counts.copy()
    for sp in counts.columns:
        own = next((g for g in LINEAGES.values() if sp in g), [sp])
        out.loc[~P.drop(columns=own).any(axis=1), sp] = 0
    return out


def lunularia_input():
    fa = ROOT / "data/orthofinder_input/Lcruciata.fa"
    L = fasta_lengths(fa)
    og = pd.read_csv(OF / "Orthogroups/Orthogroups.tsv", sep="\t", index_col=0)
    genes = {g for v in og.Lcruciata.dropna() for g in v.split(", ")}
    un = [l for g, l in L.items() if g not in genes]
    asg = [l for g, l in L.items() if g in genes]
    return pd.DataFrame([{"gene_models_in_input": len(L), "assigned_to_orthogroups": len(asg),
                          "unassigned": len(un), "fraction_unassigned": len(un) / len(L),
                          "median_length_unassigned_aa": float(np.median(un)),
                          "median_length_assigned_aa": float(np.median(asg)),
                          "published_gene_number_Linde2023": 16768}])


def sign_convention(bm, counts):
    """Pearson r between diversity and the second eigenvector, before the sign is fixed."""
    rows = []
    for label, M in (("presence", bm.build_raw_matrix(counts)), ("RCA", bm.build_rca_matrix(counts))):
        M = bm.filter_matrix(M)
        for subset in ("all eleven", "without Creinhardtii"):
            m = M if subset == "all eleven" else M.drop(columns="Creinhardtii")
            m = m[(m.sum(axis=1) >= 2)]
            mat = m.values.astype(float)
            k_sp, k_og = mat.sum(axis=0), mat.sum(axis=1)
            T = (mat.T / k_sp[:, None]) @ (mat / k_og[:, None])
            w, v = np.linalg.eig(T)
            vec = v[:, np.argsort(w.real)[::-1][1]].real
            r = np.corrcoef(k_sp, vec)[0, 1]
            vec = vec if r >= 0 else -vec
            s = pd.Series(vec, index=m.columns)
            rows.append({"matrix": label, "species": subset, "abs_r_diversity_vs_axis": abs(r),
                         "liverworts_mean": s[LIVERWORTS].mean(), "angiosperms_mean": s[ANGIOSPERMS].mean(),
                         "liverworts_and_angiosperms_at_opposite_ends":
                             bool(np.sign(s[LIVERWORTS].mean()) != np.sign(s[ANGIOSPERMS].mean()))})
    return pd.DataFrame(rows)


GO_DIR = ROOT / "CODE/data/go"


def missing_go_enrichment(counts):
    """GO biological process terms enriched among the orthogroups that both liverworts lack.

    Target: orthogroups present in all three angiosperms and in at least one other non-seed
    land plant but absent from both liverworts. Background: all orthogroups present in all
    three angiosperms and in at least one other non-seed land plant. Terms are those of the
    Arabidopsis members (EBI GOA, NOT excluded, propagated through is_a and part_of);
    hypergeometric test, Benjamini-Hochberg correction over terms with at least ten
    background orthogroups. Needs the GOA and go-basic files of CODE/data/go.
    """
    import gzip
    import sys
    from scipy.stats import hypergeom
    sys.setrecursionlimit(10000)
    par, name, space, cur = {}, {}, {}, None
    for line in open(GO_DIR / "go-basic.obo"):
        line = line.strip()
        if line == "[Term]":
            cur = None
        elif line.startswith("id: GO:"):
            cur = line[4:]; par[cur] = set()
        elif cur and line.startswith("name: "):
            name[cur] = line[6:]
        elif cur and line.startswith("namespace: "):
            space[cur] = line[11:]
        elif cur and line.startswith("is_a: "):
            par[cur].add(line[6:16])
        elif cur and line.startswith("relationship: part_of "):
            par[cur].add(line[22:32])
    anc = {}

    def ancestors(t):
        if t not in anc:
            out = {t}
            for q in par.get(t, ()):
                out |= ancestors(q)
            anc[t] = out
        return anc[t]

    gene_terms, gene_name, gene_evidence = {}, {}, {}
    for line in gzip.open(GO_DIR / "goa_arabidopsis.gaf.gz", "rt"):
        if line[0] == "!":
            continue
        f = line.rstrip("\n").split("\t")
        m = re.search(r"\bAT[1-5MC]G\d{5}\b", f[10], re.I)       # TAIR locus lives in column 11
        if not m:
            continue
        g = m.group(0).upper()
        gene_name[g] = f[9]
        if "NOT" not in f[3]:
            gene_terms.setdefault(g, set()).add(f[4])
            gene_evidence.setdefault(g, []).append((f[4], f[6]))
    members = pd.read_csv(OF / "Orthogroups/Orthogroups.tsv", sep="\t", index_col=0).Athaliana
    og_terms, og_ppr = {}, {}
    for o in counts.index:
        v = members.get(o)
        if not isinstance(v, str):
            continue
        genes = [x.split(".")[0].upper() for x in v.split(", ")]
        og_terms[o] = {a for g in genes for t in gene_terms.get(g, ()) if t in par for a in ancestors(t)}
        og_ppr[o] = any(re.search("pentatricopeptide", gene_name.get(g, ""), re.I) for g in genes)
    P = counts > 0
    conserved = P[ANGIOSPERMS].all(axis=1) & P[OTHER_NON_SEED].any(axis=1)
    missing = conserved & ~P[LIVERWORTS].any(axis=1)
    bg = [o for o in counts.index[conserved] if o in og_terms]
    tg = [o for o in counts.index[missing] if o in og_terms]
    from collections import Counter
    cb = Counter(t for o in bg for t in og_terms[o])
    ct = Counter(t for o in tg for t in og_terms[o])
    rows = []
    for t, K in cb.items():                       # every term with at least ten background orthogroups
        k = ct.get(t, 0)
        if K < 10 or space.get(t) != "biological_process":
            continue
        hits = [o for o in tg if t in og_terms[o]]
        rows.append({"GO": t, "term": name[t], "missing_orthogroups": k, "background_orthogroups": K,
                     "fold": k / (len(tg) * K / len(bg)), "p": hypergeom.sf(k - 1, len(bg), K, len(tg)),
                     "with_PPR_protein": sum(og_ppr[o] for o in hits)})
    d = pd.DataFrame(rows).sort_values(["p", "GO"], kind="mergesort").reset_index(drop=True)   # ties in GO order
    m = len(d)
    d["q"] = (d.p * m / np.arange(1, m + 1))[::-1].cummin()[::-1].clip(upper=1)
    d.attrs["n"] = (len(tg), len(bg))
    top = d.GO.iloc[0]
    desc = {t for t in par if top in ancestors(t)}
    rows = []
    for o in [o for o in tg if top in og_terms[o]]:
        genes = [x.split(".")[0].upper() for x in members[o].split(", ")]
        ev = [e for g in genes for t, e in gene_evidence.get(g, []) if t in desc]
        row = {"orthogroup": o, "PPR_protein": og_ppr[o],
               "evidence_IEA_or_IBA": sum(e in ("IEA", "IBA") for e in ev), "evidence_other": sum(e not in ("IEA", "IBA") for e in ev)}
        row.update({sp: int(counts.at[o, sp]) for sp in counts.columns})
        rows.append(row)
    d.attrs["details"] = pd.DataFrame(rows)
    return d


def hog_lineages(counts):
    """Liverwort genes per ancestral gene lineage of the angiosperms (node N8) and of the seed plants (N7),
    in land plant hierarchical orthogroups (N1) held by all three angiosperms and both liverworts,
    separately for orthogroups with and without an Arabidopsis transcription factor (PlantTFDB)."""
    H = TREES / "Phylogenetic_Hierarchical_Orthogroups"
    tf = set(pd.read_csv(TF_LIST, sep="\t").Gene_ID)
    genes = lambda v: [x.strip() for x in v.split(",")] if isinstance(v, str) and v.strip() else []
    n1 = pd.read_csv(H / "N1.tsv", sep="\t", dtype=str)
    sp = [c for c in n1.columns if c not in ("HOG", "OG", "Gene Tree Parent Clade")]
    rec, g2n1 = {}, {}
    for _, r in n1.iterrows():
        rec[r.HOG] = {x: genes(r[x]) for x in sp}
        for x in sp:
            for g in rec[r.HOG][x]:
                g2n1[g] = r.HOG

    def lineages(node, members):
        c = {}
        for _, r in pd.read_csv(H / f"{node}.tsv", sep="\t", dtype=str).iterrows():
            for h in {g2n1.get(g) for x in members for g in genes(r[x])} - {None}:
                c[h] = c.get(h, 0) + 1
        return c

    c8, c6 = lineages("N8", ANGIOSPERMS), lineages("N7", ANGIOSPERMS + ["Pabies"])
    rows = []
    for h, o in rec.items():
        if all(o[a] for a in ANGIOSPERMS) and all(o[l] for l in LIVERWORTS):
            rows.append({"tf": any(g.split(".")[0] in tf for g in o["Athaliana"]),
                         "liverwort_genes": np.mean([len(o[l]) for l in LIVERWORTS]),
                         "angiosperm_genes": np.mean([len(o[a]) for a in ANGIOSPERMS]),
                         "N8": c8.get(h, 0), "N7": c6.get(h, 0)})
    df = pd.DataFrame(rows)
    out = []
    for k, g in df.groupby("tf"):
        out.append({"orthogroups": "with a transcription factor" if k else "without a transcription factor",
                    "hierarchical_orthogroups": len(g),
                    "liverwort_genes_per_angiosperm_lineage_N8": g.liverwort_genes.sum() / g.N8.sum(),
                    "liverwort_genes_per_seed_plant_lineage_N7": g.liverwort_genes.sum() / g.N7.sum(),
                    "angiosperm_genes_per_angiosperm_lineage_N8": g.angiosperm_genes.sum() / g.N8.sum()})
    return pd.DataFrame(out)


def dyw_counts():
    """Proteins ending in the canonical DYW domain (C-x-C-x(1,2)-D-Y-W in the last 12 residues) per species,
    and the orthogroups that hold them, with the number of angiosperm genes in each."""
    og = pd.read_csv(OF / "Orthogroups/Orthogroups.tsv", sep="\t", index_col=0)
    gene2og = {g.strip(): o for o, r in og.iterrows() for s in og.columns if isinstance(r[s], str)
               for g in r[s].split(", ")}
    rows, ogs = [], {}
    for fa in sorted((ROOT / "data/orthofinder_input").glob("*.fa")):
        sp = fa.stem
        L = fasta_lengths(fa)                         # names normalised as in Orthogroups.tsv
        seqs, k = {}, None
        for line in open(fa):
            if line.startswith(">"):
                k = line[1:].split()[0].replace(":", "_")
                seqs[k] = []
            else:
                seqs[k].append(line.strip().rstrip("*"))
        hits = [g for g, v in seqs.items() if re.search(r"C.C..?DYW$", "".join(v)[-12:])]
        rows.append({"species": sp, "DYW_proteins": len(hits)})
        for g in hits:
            o = gene2og.get(g)
            if o:
                ogs.setdefault(o, {}).setdefault(sp, 0)
                ogs[o][sp] += 1
    per_og = []
    for o, d in ogs.items():
        r = {"orthogroup": o, **d}
        r["angiosperm_genes_in_orthogroup"] = sum(len(og.at[o, a].split(", ")) for a in ANGIOSPERMS
                                                  if isinstance(og.at[o, a], str))
        per_og.append(r)
    return pd.DataFrame(rows), pd.DataFrame(per_og).fillna(0)


TF_LIST = ROOT / "data/tf_families/Ath_TF_list.txt"


STRATA = [(0, 1, "up to one"), (1, 2, "one to two"), (2, 4, "two to four"), (4, np.inf, "more than four")]
PAIRS = [("KEGG regulatory", "KEGG primary"), ("GO regulatory", "GO primary"),
         ("transcription factors", "KEGG primary")]


def regulatory_stratified(counts):
    """Median liverwort/angiosperm copy ratio of regulatory and metabolic orthogroups within strata of
    mean angiosperm copy number, with the orthogroups shared by both categories removed (as in 27)."""
    cat = pd.read_csv(ROOT / "data/results_regulatory_contrast/category_orthogroups.tsv", sep="\t")
    sets = {k: set(v) for k, v in cat.groupby("category").orthogroup}
    occ_mean = lambda x: x.where(x > 0).mean(axis=1)
    ok = (counts[LIVERWORTS] > 0).any(axis=1) & (counts[ANGIOSPERMS] > 0).any(axis=1)
    ratio = (occ_mean(counts[LIVERWORTS]) / occ_mean(counts[ANGIOSPERMS]))[ok]
    angio = occ_mean(counts[ANGIOSPERMS])[ok]
    from scipy.stats import mannwhitneyu
    rows = []
    for reg, met in PAIRS:
        members = {reg: sets[reg] - sets[met], met: sets[met] - sets[reg]}
        idx = {k: ratio.index.intersection(list(v)) for k, v in members.items()}
        stratum = {k: {name: idx[k][(angio[idx[k]] > lo) & (angio[idx[k]] <= hi)] for lo, hi, name in STRATA}
                   for k in idx}
        for label in (reg, met):
            row = {"comparison": f"{reg} vs {met}", "category": label, "orthogroups": len(idx[label]),
                   "median_angiosperm_copies": angio[idx[label]].median()}
            for lo, hi, name in STRATA:
                sel = stratum[label][name]
                row[f"n {name}"] = len(sel)
                row[f"median ratio {name}"] = ratio[sel].median() if len(sel) else np.nan
                a, b = ratio[stratum[reg][name]], ratio[stratum[met][name]]
                # two-sided Mann-Whitney test of the regulatory against the metabolic ratios in the stratum
                row[f"P {name}"] = (mannwhitneyu(a, b, alternative="two-sided").pvalue
                                    if label == reg and len(a) and len(b) else np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    bm = load_bm()
    all_counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    counts = all_counts[(all_counts > 0).sum(axis=1) >= 2]  # same min_species filter as 05
    if "--sign-only" in __import__("sys").argv:
        sign_convention(bm, all_counts).to_csv(OUT / "M_sign_convention.tsv", sep="\t", index=False, float_format="%.4f")
        return
    if "--stratified-only" in __import__("sys").argv:
        strat = regulatory_stratified(counts)
        strat.to_csv(OUT / "M_regulatory_stratified.tsv", sep="\t", index=False, float_format="%.4g")
        print(strat.round(3).to_string(index=False))
        return

    ov = overlap(counts)
    base = filterB_baseline()

    sens = []
    for label, c in (("orthogroups (reported)", counts), ("hierarchical orthogroups, land plant node (N1)",
                                                          hog_counts("N1"))):
        sens.append((label, *ratios(c)))
    frag, _ = without_fragments(bm)
    sens.append(("genes < 0.5 x orthogroup median length removed", *ratios(frag)))
    sens.append(("lineage-specific orthogroups removed", *ratios(without_lineage_specific(counts))))
    sens = pd.DataFrame(sens, columns=["set", "coverage_ratio_L_A", "genes_per_OG_ratio_L_A"])

    lun = lunularia_input()
    sign = sign_convention(bm, all_counts)                 # RCA over all orthogroups, then filtered, as in 05
    go = missing_go_enrichment(counts) if (GO_DIR / "go-basic.obo").exists() else None
    losses = lineage_losses(counts)
    hogs = hog_lineages(counts) if TF_LIST.exists() else None
    dyw, dyw_og = dyw_counts()
    deficit = deficit_decomposition(counts)
    expected = busco_expected(ov, losses)
    best = kegg_best_hit()

    OUT.mkdir(parents=True, exist_ok=True)
    ov.to_csv(OUT / "M_overlap_with_angiosperms.tsv", sep="\t", index=False, float_format="%.6f")
    base.to_csv(OUT / "M_filterB_sharing_baseline.tsv", sep="\t", index=False, float_format="%.4f")
    sens.to_csv(OUT / "M_sensitivity_ratios.tsv", sep="\t", index=False, float_format="%.4f")
    lun.to_csv(OUT / "M_lunularia_input.tsv", sep="\t", index=False, float_format="%.4f")
    sign.to_csv(OUT / "M_sign_convention.tsv", sep="\t", index=False, float_format="%.4f")
    losses.to_csv(OUT / "M_lineage_losses.tsv", sep="\t", index=False)
    dyw.to_csv(OUT / "M_dyw_counts.tsv", sep="\t", index=False)
    dyw_og.to_csv(OUT / "M_dyw_orthogroups.tsv", sep="\t", index=False)
    print("\nDYW proteins\n" + dyw.to_string(index=False))
    print(dyw_og[dyw_og.get("Ppatens", 0) > 0].to_string(index=False) if "Ppatens" in dyw_og else "")
    if hogs is not None:
        hogs.to_csv(OUT / "M_hog_lineages.tsv", sep="\t", index=False, float_format="%.4f")
        print("\nhierarchical orthogroup lineages\n" + hogs.round(3).to_string(index=False))
    expected.to_csv(OUT / "M_busco_expected_loss.tsv", sep="\t", index=False, float_format="%.4f")
    print("\nloss expected from BUSCO\n" + expected[expected.species.isin(LIVERWORTS + ANGIOSPERMS)].round(4).to_string(index=False))
    if best is not None:
        best.to_csv(OUT / "M_kegg_best_hit.tsv", sep="\t", index=False)
        print("\nbest-hit KO sensitivity\n" + best.to_string(index=False))
    deficit.to_csv(OUT / "M_deficit_decomposition.tsv", sep="\t", index=False, float_format="%.4f")
    regulatory_stratified(counts).to_csv(OUT / "M_regulatory_stratified.tsv", sep="\t", index=False,
                                         float_format="%.4g")
    print("\nlineage losses\n" + losses.to_string(index=False))
    print("\ndeficit decomposition\n" + deficit.round(3).to_string(index=False))
    if go is not None:
        go.to_csv(OUT / "M_missing_GO_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
        go.attrs["details"].to_csv(OUT / "M_missing_GO_top_term_orthogroups.tsv", sep="\t", index=False)
        print(f"\nGO enrichment of orthogroups missing from both liverworts (n = {go.attrs['n'][0]} of "
              f"{go.attrs['n'][1]})\n" + go.head(8).round(4).to_string(index=False))

    pd.set_option("display.width", 220)
    for name, d in (("overlap", ov), ("Filter B baseline", base), ("sensitivity", sens),
                    ("Lunularia input", lun), ("sign convention", sign)):
        print(f"\n{name}\n{d.round(3).to_string(index=False)}")


if __name__ == "__main__":
    main()
