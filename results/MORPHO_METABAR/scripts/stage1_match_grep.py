#!/usr/bin/env python3
"""
Stage 1 — Match-grep morfología ↔ linajes BOLD (EDA_TAX_RUN29/tax_parsed.tsv)

Estrategia:
  a) Match EXACTO por token contra cada rango taxonómico (Phylum, Class, Order,
     Family, Genus) del tax_parsed BOLD.
  b) Match RELAJADO por substring (`str.contains`, case-insensitive) para
     capturar variantes (Siphonophorae ↔ Siphonophora, Hydromedusae ↔ Hydrozoa,
     Larva_briozoa ↔ Bryozoa, etc.).

Salidas:
  reports/stage1_match_report.tsv
  reports/stage1_summary.md
"""
import re
import sys
from pathlib import Path
import pandas as pd

LANCO = Path("/sessions/gracious-sleepy-edison/mnt/LANCO")
TAX   = LANCO / "results/EDA_TAX_RUN29/tables/tax_parsed.tsv"
MORPHO_LONG = LANCO / "results/MORPHO_METABAR/tables/morpho_long.tsv"
OUT = LANCO / "results/MORPHO_METABAR/reports"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Lista morfológica (tal como fue reportada por el laboratorio)
# ---------------------------------------------------------------------------
MORPHO_TAXA = [
    "Calanoida", "Cyclopoida", "Poecilostomatoida", "Cladocera", "Mysidacea",
    "Nauplio", "Zoea", "Megalopa", "Euphausiidae_adulto", "Larva_briozoa",
    "Larva_veliger_Atlantidae", "Larva_veliger", "Larva_cirripeda",
    "Larva_bivalvia", "Larva_polychaeta", "Larva_pez", "Huevos_pez",
    "Chaetognatha", "Siphonophorae", "Hydromedusae", "Doliolida",
    "Appendicularia", "Nematoda", "Platyhelminthes",
]

# ---------------------------------------------------------------------------
# Cargar linaje BOLD parseado (261 ASVs × 8 rangos)
# ---------------------------------------------------------------------------
tax = pd.read_csv(TAX, sep="\t")
RANKS = ["Phylum", "Class", "Order", "Family", "Genus"]

# Sets únicos por rango (para grep rápido)
rank_sets = {r: set(tax[r].dropna().unique()) for r in RANKS}
rank_str  = {r: "|".join(sorted(rank_sets[r])) for r in RANKS}

# ---------------------------------------------------------------------------
# Núcleo del match: normalización + búsqueda
# ---------------------------------------------------------------------------
def _clean(name: str) -> str:
    """Extraer palabra biológica útil del nombre morfológico"""
    # Quitar prefijos larvarios / sufijos de estadio
    n = re.sub(r"^(Larva_|Huevos_)", "", name, flags=re.I)
    n = re.sub(r"_(adulto|larva|juvenil)$", "", n, flags=re.I)
    n = n.replace("_", " ").strip()
    return n

def match_exact(query: str):
    """Coincidencia estricta case-insensitive por rango"""
    q = query.lower()
    hits = {}
    for r in RANKS:
        matched = [t for t in rank_sets[r] if isinstance(t, str) and t.lower() == q]
        if matched:
            hits[r] = matched
    return hits

def match_relaxed(query: str):
    """Contains: query in taxon o taxon in query, case-insensitive"""
    q = query.lower()
    hits = {}
    for r in RANKS:
        matched = [t for t in rank_sets[r] if isinstance(t, str)
                   and (q in t.lower() or t.lower() in q)]
        if matched:
            hits[r] = matched
    return hits

# ---------------------------------------------------------------------------
# Ejecutar match y consolidar
# ---------------------------------------------------------------------------
rows = []
for name in MORPHO_TAXA:
    cleaned = _clean(name)
    ex = match_exact(cleaned)
    rl = match_relaxed(cleaned)
    # Nivel más informativo del match exacto
    if ex:
        rank_hit_exact = next((r for r in RANKS if r in ex), None)
        taxa_exact = ";".join(ex.get(rank_hit_exact, []))
    else:
        rank_hit_exact, taxa_exact = None, ""
    # Nivel más informativo del match relajado
    if rl:
        rank_hit_relax = next((r for r in RANKS if r in rl), None)
        taxa_relax = ";".join(sorted(set(sum(rl.values(), []))))
    else:
        rank_hit_relax, taxa_relax = None, ""

    # Contar ASVs BOLD asociados al match relajado (por cualquier rango)
    if rl:
        mask = pd.Series(False, index=tax.index)
        for r, ts in rl.items():
            mask |= tax[r].isin(ts)
        n_asvs = int(mask.sum())
        asv_ids = ";".join(tax.loc[mask, "ASV_ID"].head(5).tolist())
    else:
        n_asvs, asv_ids = 0, ""
    rows.append({
        "morpho_taxon": name,
        "morpho_normalized": cleaned,
        "match_exact_rank": rank_hit_exact,
        "match_exact_taxa": taxa_exact,
        "match_relaxed_rank": rank_hit_relax,
        "match_relaxed_taxa": taxa_relax,
        "n_asvs_bold_hit": n_asvs,
        "example_asvs": asv_ids,
        "status": "EXACT" if ex else ("RELAXED" if rl else "NO_MATCH"),
    })

report = pd.DataFrame(rows)
report.to_csv(OUT / "stage1_match_report.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Resumen textual
# ---------------------------------------------------------------------------
counts = report.status.value_counts().to_dict()
lines = [
    "# Stage 1 — match-grep morfología ↔ BOLD",
    "",
    f"- Taxa morfológicos evaluados: **{len(MORPHO_TAXA)}**",
    f"- Match exacto: **{counts.get('EXACT', 0)}**",
    f"- Match relajado (substring): **{counts.get('RELAXED', 0)}**",
    f"- Sin coincidencia: **{counts.get('NO_MATCH', 0)}**",
    "",
    "## Tabla-resumen",
    "",
    "| Morfológico | Normalizado | Estado | Rango BOLD | Taxa BOLD hit | ASVs |",
    "|---|---|---|---|---|---|",
]
for _, r in report.iterrows():
    rank = r["match_exact_rank"] or r["match_relaxed_rank"] or "—"
    taxa = r["match_exact_taxa"] or r["match_relaxed_taxa"] or "—"
    lines.append(
        f"| {r['morpho_taxon']} | {r['morpho_normalized']} | {r['status']} | "
        f"{rank or '—'} | {taxa or '—'} | {r['n_asvs_bold_hit']} |"
    )

(OUT / "stage1_summary.md").write_text("\n".join(lines))
print("\n".join(lines[:10]))
print(f"\nSalidas escritas en: {OUT}")
print(report[['morpho_taxon','status','match_exact_rank','match_relaxed_rank','n_asvs_bold_hit']].to_string(index=False))
