#!/usr/bin/env python3
"""
Stage 2 — Búsqueda a fondo en el clasificador BOLD (todos los rangos, incluida
especie y BOLD BIN) y en los hits BLAST completos.

Fuentes locales inspeccionadas:
  - results/EDA_TAX_RUN29/tables/tax_parsed.tsv       (261 ASVs × 8 rangos)
  - results/04_taxonomy/tax_BOLD.tsv                  (bruto BOLD)
  - results/04_taxonomy/tax_BOLD_search.qza           (blast6.tsv con todos los hits)
  - results/04_taxonomy/tax_BOLD_classification.qza   (taxonomy.tsv oficial)

Bases de referencia globales (`db/rescript/bold-coi-tax.qza` y
`db/rescript/ncbi-coi-tax.qza`) NO están en el repo local — están en LUSTRE.
Ver bloque final para el comando esperado si se ejecuta en el cluster.
"""
import io
import re
import zipfile
from pathlib import Path
import pandas as pd

LANCO = Path("/sessions/gracious-sleepy-edison/mnt/LANCO")
TAX_PARSED = LANCO / "results/EDA_TAX_RUN29/tables/tax_parsed.tsv"
TAX_BOLD   = LANCO / "results/04_taxonomy/tax_BOLD.tsv"
SEARCH_QZA = LANCO / "results/04_taxonomy/tax_BOLD_search.qza"
OUT = LANCO / "results/MORPHO_METABAR/reports"
OUT.mkdir(parents=True, exist_ok=True)

MORPHO_TAXA = [
    "Calanoida", "Cyclopoida", "Poecilostomatoida", "Cladocera", "Mysidacea",
    "Nauplio", "Zoea", "Megalopa", "Euphausiidae_adulto", "Larva_briozoa",
    "Larva_veliger_Atlantidae", "Larva_veliger", "Larva_cirripeda",
    "Larva_bivalvia", "Larva_polychaeta", "Larva_pez", "Huevos_pez",
    "Chaetognatha", "Siphonophorae", "Hydromedusae", "Doliolida",
    "Appendicularia", "Nematoda", "Platyhelminthes",
]

# Variantes de root-word (a partir del *nombre morfológico* -> raíces biológicas
# esperables en linajes BOLD/NCBI, sin mapear todavía la nomenclatura completa
# WoRMS). Solo raíces léxicas — no es aún el mapping WoRMS de Stage 3.
LEXICAL_ROOTS = {
    "Calanoida":               ["calanoid"],
    "Cyclopoida":              ["cyclopoid"],
    "Poecilostomatoida":       ["poecilostomat", "poecilo"],
    "Cladocera":               ["cladocer", "branchiopod"],
    "Mysidacea":               ["mysid"],
    "Nauplio":                 ["copepod", "cirripedia", "nauplius"],
    "Zoea":                    ["decapod", "brachyur", "zoea"],
    "Megalopa":                ["brachyur", "decapod", "megalop"],
    "Euphausiidae_adulto":     ["euphaus"],
    "Larva_briozoa":           ["bryozoa", "briozoa", "gymnolaem", "cheilostomat"],
    "Larva_veliger_Atlantidae":["atlantidae", "atlanta", "carinari", "pterotrach"],
    "Larva_veliger":           ["mollusc", "veliger", "gastropod", "bivalv"],
    "Larva_cirripeda":         ["cirripedia", "balan", "thoracica"],
    "Larva_bivalvia":          ["bivalv", "pelecypod", "ostrea", "mytil"],
    "Larva_polychaeta":        ["polychaet", "annelid"],
    "Larva_pez":               ["teleost", "actinopter", "chordat", "fish"],
    "Huevos_pez":              ["teleost", "actinopter"],
    "Chaetognatha":            ["chaetognath", "sagittoid", "sagitta"],
    "Siphonophorae":           ["siphonoph"],
    "Hydromedusae":            ["hydroz", "hydromedus", "leptothec", "anthomedus"],
    "Doliolida":               ["doliolid", "thaliace"],
    "Appendicularia":          ["appendicul", "oikopleur", "larvace"],
    "Nematoda":                ["nematod"],
    "Platyhelminthes":         ["platyhelminth", "turbellar", "trematod", "cestod"],
}

# ──────────────────────────────────────────────────────────────────────────
# Cargar linajes conocidos localmente
# ──────────────────────────────────────────────────────────────────────────
tax_parsed = pd.read_csv(TAX_PARSED, sep="\t")
RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "BOLD_BIN"]

# Cadena concatenada por ASV para búsqueda substring en TODOS los rangos
tax_parsed["_concat"] = tax_parsed[RANKS].fillna("").agg(";".join, axis=1)

# ──────────────────────────────────────────────────────────────────────────
# Cargar blast6.tsv completo desde el QZA de búsqueda (todos los hits)
# ──────────────────────────────────────────────────────────────────────────
with zipfile.ZipFile(SEARCH_QZA) as z:
    hit_name = [n for n in z.namelist() if n.endswith("data/blast6.tsv")][0]
    with z.open(hit_name) as fh:
        blast = pd.read_csv(fh, sep="\t", header=None,
            names=["qseqid","sseqid","pident","length","mismatch","gapopen",
                   "qstart","qend","sstart","send","evalue","bitscore"])
n_hits = len(blast)
n_asvs_hit = blast["qseqid"].nunique()
n_bold_refs = blast["sseqid"].nunique()

# ──────────────────────────────────────────────────────────────────────────
# Buscar cada taxón morfológico por raíces léxicas en linajes locales
# ──────────────────────────────────────────────────────────────────────────
def scan_local(morpho_name):
    roots = LEXICAL_ROOTS[morpho_name]
    concat_lc = tax_parsed["_concat"].str.lower()
    hits = pd.Series(False, index=tax_parsed.index)
    matched_roots = []
    for root in roots:
        m = concat_lc.str.contains(re.escape(root.lower()), na=False)
        if m.any():
            matched_roots.append(root)
            hits |= m
    n = int(hits.sum())
    example_lineage = tax_parsed.loc[hits, "_concat"].head(3).tolist()
    return {
        "morpho_taxon": morpho_name,
        "lexical_roots_tried": ";".join(roots),
        "roots_found": ";".join(matched_roots) if matched_roots else "",
        "n_asvs_with_root": n,
        "example_lineages": " || ".join(example_lineage),
        "present_in_tax_parsed": n > 0,
    }

report = pd.DataFrame([scan_local(m) for m in MORPHO_TAXA])
report.to_csv(OUT / "stage2_db_scan_report.tsv", sep="\t", index=False)

# ──────────────────────────────────────────────────────────────────────────
# Cruzar con Stage 1 para consolidar estado
# ──────────────────────────────────────────────────────────────────────────
stage1 = pd.read_csv(OUT / "stage1_match_report.tsv", sep="\t")
merge  = report.merge(stage1[["morpho_taxon", "status"]], on="morpho_taxon")
merge  = merge.rename(columns={"status": "stage1_status"})
merge["stage2_status"] = merge["n_asvs_with_root"].apply(
    lambda n: "PRESENT_LOCAL" if n > 0 else "ABSENT_LOCAL")
merge.to_csv(OUT / "stage2_consolidated.tsv", sep="\t", index=False)

# ──────────────────────────────────────────────────────────────────────────
# Reporte legible
# ──────────────────────────────────────────────────────────────────────────
lines = [
    "# Stage 2 — Presencia de morfología en linajes BOLD locales",
    "",
    f"- Registros en `tax_BOLD_search.qza`/blast6: **{n_hits}** hits sobre "
    f"**{n_asvs_hit}** ASVs y **{n_bold_refs}** referencias BOLD únicas",
    f"- Búsqueda léxica sobre linajes concatenados de `tax_parsed.tsv` "
    f"(8 rangos, 261 ASVs)",
    "",
    "## Resultado por taxón morfológico",
    "",
    "| Morfológico | Raíces léxicas probadas | Encontrado en linaje BOLD | ASVs | Stage 1 | Stage 2 | Ejemplo linaje |",
    "|---|---|---|---|---|---|---|",
]
for _, r in merge.iterrows():
    lines.append(
        f"| {r['morpho_taxon']} | `{r['lexical_roots_tried']}` | "
        f"`{r['roots_found'] or '—'}` | {r['n_asvs_with_root']} | "
        f"{r['stage1_status']} | {r['stage2_status']} | "
        f"{r['example_lineages'][:110] + '…' if len(r['example_lineages'])>110 else r['example_lineages']} |"
    )
lines += [
    "",
    "## Ejecución complementaria en el cluster LUSTRE (recomendada)",
    "",
    "Las bases de referencia completas están fuera del repo local. Correr:",
    "",
    "```bash",
    "# Extraer taxonomía cruda de las DBs de referencia y grepear la lista",
    "cd $LANCO/db/rescript",
    "for f in bold-coi-tax.qza ncbi-coi-tax.qza; do",
    "  unzip -p $f '*/data/taxonomy.tsv' > /tmp/${f%.qza}.tsv",
    "done",
    "# Ejemplo: grep de todas las raíces léxicas producidas por este script",
    "grep -iE 'calanoid|cyclopoid|cladocer|mysid|euphaus|bryozoa|cirripedia|polychaet|siphonoph|hydromedus|doliolid|appendicul|nematod|platyhelminth' /tmp/bold-coi-tax.tsv | wc -l",
    "```",
]
(OUT / "stage2_summary.md").write_text("\n".join(lines))
print("\n".join(lines[:10]))
print()
print(merge[["morpho_taxon","stage1_status","stage2_status","n_asvs_with_root","roots_found"]].to_string(index=False))
