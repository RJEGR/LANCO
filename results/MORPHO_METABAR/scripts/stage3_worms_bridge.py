#!/usr/bin/env python3
"""
Stage 3 — Puente WoRMS ↔ BOLD/NCBI para reconciliar nomenclatura morfológica.

Estrategia en 3 capas:
  A) MAPPING CURADO (offline). Cada nombre morfológico se resuelve a su
     equivalente taxonómico Linneano en WoRMS/NCBI, incluyendo:
       - taxón_referencia_worms  (nombre aceptado en WoRMS)
       - rango_referencia        (rango en el que ese nombre existe en BOLD)
       - notas                   (razón del mapping: subsunción, larva, sinonimia)

  B) VERIFICACIÓN CONTRA LINAJES BOLD LOCAL (tax_parsed.tsv). Para cada
     mapping, se comprueba si el taxón_referencia aparece en algún rango de la
     tabla de linajes.

  C) API WoRMS (opcional). Si hay red, resuelve cada nombre contra
     https://www.marinespecies.org/rest y obtiene el AphiaID + classification
     canónica. Fallback silencioso a la capa A si la red no está disponible.

Salidas:
  reports/stage3_worms_bridge.tsv           mapping curado + verificación
  reports/stage3_final_convergence.tsv      tabla final morfología ↔ BOLD
  reports/stage3_summary.md                 informe legible
"""
import json
import re
import urllib.request
import urllib.parse
import urllib.error
from pathlib import Path
import pandas as pd

LANCO = Path("/sessions/gracious-sleepy-edison/mnt/LANCO")
TAX_PARSED = LANCO / "results/EDA_TAX_RUN29/tables/tax_parsed.tsv"
OUT = LANCO / "results/MORPHO_METABAR/reports"

# ──────────────────────────────────────────────────────────────────────────
# A) MAPPING CURADO (basado en WoRMS + Boxshall & Halsey 2004 + literatura)
# ──────────────────────────────────────────────────────────────────────────
# Formato: morpho -> [(worms_name, target_rank_in_BOLD, note)]
# Si un nombre morfológico agrupa >1 taxón (p. ej. Nauplio = larva de Copepoda
# ∪ Cirripedia), se listan todos los targets.
CURATED = {
    "Calanoida":                [("Calanoida",       "Order",  "Aceptado WoRMS")],
    "Cyclopoida":               [("Cyclopoida",      "Order",  "Aceptado WoRMS")],
    "Poecilostomatoida":        [("Cyclopoida",      "Order",  "Boxshall & Halsey 2004: Poecilostomatoida subsumido en Cyclopoida")],
    "Cladocera":                [("Branchiopoda",    "Class",  "Cladocera es superorden dentro de Branchiopoda; BOLD usa Class Branchiopoda")],
    "Mysidacea":                [("Mysida",          "Order",  "Mysidacea antiguo; WoRMS acepta Order Mysida")],
    "Nauplio":                  [("Copepoda",        "Class",  "Estadio larvario copépodo"),
                                 ("Thecostraca",     "Class",  "Estadio larvario cirripedios (Thecostraca)")],
    "Zoea":                     [("Decapoda",        "Order",  "Estadio zoea de Decapoda"),
                                 ("Malacostraca",    "Class",  "Ancestro común de Decapoda")],
    "Megalopa":                 [("Brachyura",       "Infraorder", "Megalopa = estadio postlarva de Brachyura (dentro de Decapoda)"),
                                 ("Decapoda",        "Order",  "BOLD colapsa a Order Decapoda")],
    "Euphausiidae_adulto":      [("Euphausiacea",    "Order",  "Familia Euphausiidae → único Order Euphausiacea"),
                                 ("Euphausiidae",    "Family", "Nombre morfológico ya es Family")],
    "Larva_briozoa":            [("Bryozoa",         "Phylum", "Cifonauta/coronado, larvas planctónicas de Bryozoa")],
    "Larva_veliger_Atlantidae": [("Atlantidae",      "Family", "Veliger heteropodo de Atlantidae (Gastropoda pelágica)")],
    "Larva_veliger":            [("Gastropoda",      "Class",  "Veliger genérico; sin diferenciar Gastropoda vs Bivalvia"),
                                 ("Bivalvia",        "Class",  "Ídem si es bivalvo")],
    "Larva_cirripeda":          [("Thecostraca",     "Class",  "Nauplius/cyprid de Cirripedia (Thecostraca)"),
                                 ("Cirripedia",      "Infraclass", "Nombre clásico aún usado")],
    "Larva_bivalvia":           [("Bivalvia",        "Class",  "Veliger de bivalvos")],
    "Larva_polychaeta":         [("Polychaeta",      "Class",  "Trocófora / nectóquete de poliquetos")],
    "Larva_pez":                [("Actinopterygii",  "Class",  "Larvas y juveniles de peces óseos"),
                                 ("Teleostei",       "Class",  "BOLD marca Teleostei como Class")],
    "Huevos_pez":               [("Actinopterygii",  "Class",  "Huevos flotantes de teleósteos"),
                                 ("Teleostei",       "Class",  "Ídem")],
    "Chaetognatha":             [("Chaetognatha",    "Phylum", "Aceptado WoRMS"),
                                 ("Sagittoidea",     "Class",  "BOLD usa Class Sagittoidea")],
    "Siphonophorae":            [("Siphonophorae",   "Order",  "Aceptado WoRMS (Cnidaria: Hydrozoa)")],
    "Hydromedusae":             [("Hydrozoa",        "Class",  "Nombre morfológico para medusas de Hydrozoa")],
    "Doliolida":                [("Doliolida",       "Order",  "Aceptado WoRMS (Tunicata: Thaliacea)"),
                                 ("Thaliacea",       "Class",  "BOLD puede colapsar a Class Thaliacea")],
    "Appendicularia":           [("Appendicularia",  "Class",  "= Larvacea; Tunicata")],
    "Nematoda":                 [("Nematoda",        "Phylum", "Aceptado WoRMS")],
    "Platyhelminthes":          [("Platyhelminthes", "Phylum", "Aceptado WoRMS")],
}

# ──────────────────────────────────────────────────────────────────────────
# B) VERIFICACIÓN CONTRA LINAJES BOLD LOCAL
# ──────────────────────────────────────────────────────────────────────────
tax = pd.read_csv(TAX_PARSED, sep="\t")
RANKS = ["Kingdom","Phylum","Class","Order","Family","Genus","Species","BOLD_BIN"]

tax_concat_series = tax[RANKS].fillna("").apply(lambda r: ";".join(r.astype(str)), axis=1)
tax_concat_lc = tax_concat_series.str.lower()

def verify_in_bold(name):
    """Devuelve (nº ASVs con `name` en cualquier rango, ejemplo linaje)"""
    m = tax_concat_lc.str.contains(rf"\b{re.escape(name.lower())}\b", regex=True, na=False)
    n = int(m.sum())
    ex = tax_concat_series.loc[m].head(2).tolist()
    return n, " || ".join(ex)

# ──────────────────────────────────────────────────────────────────────────
# C) WoRMS REST API (opcional, con timeout defensivo)
# ──────────────────────────────────────────────────────────────────────────
def worms_lookup(name, timeout=8):
    """Intenta resolver `name` en WoRMS. Devuelve dict con AphiaID + classification."""
    try:
        url = (f"https://www.marinespecies.org/rest/AphiaIDByName/"
               f"{urllib.parse.quote(name)}?marine_only=false")
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            aphia = resp.read().decode().strip()
            if not aphia or aphia in ("-999", "null"):
                return None
        aphia = aphia.strip('"')
        url = f"https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/{aphia}"
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            classif = json.loads(resp.read().decode())
        # Aplanar classification
        flat = []
        cur = classif
        while cur:
            flat.append({"rank": cur.get("rank"), "scientificname": cur.get("scientificname")})
            cur = cur.get("child")
        return {"AphiaID": aphia, "classification": flat}
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError,
            json.JSONDecodeError, ConnectionError, OSError):
        return None

# ──────────────────────────────────────────────────────────────────────────
# Ejecutar mapping + verificación
# ──────────────────────────────────────────────────────────────────────────
rows = []
for morpho, mappings in CURATED.items():
    for worms_name, target_rank, note in mappings:
        n_asvs, example = verify_in_bold(worms_name)
        worms_hit = worms_lookup(worms_name)  # None si no hay red
        rows.append({
            "morpho_taxon":     morpho,
            "worms_name":       worms_name,
            "target_rank_bold": target_rank,
            "note":             note,
            "verified_in_bold": n_asvs > 0,
            "n_asvs_bold":      n_asvs,
            "example_lineage":  example[:150],
            "worms_aphia_id":   (worms_hit or {}).get("AphiaID", ""),
            "worms_classification": (
                ";".join(f"{c['rank']}:{c['scientificname']}"
                         for c in (worms_hit or {}).get("classification", []))
                if worms_hit else ""),
        })

bridge = pd.DataFrame(rows)
bridge.to_csv(OUT / "stage3_worms_bridge.tsv", sep="\t", index=False)

# ──────────────────────────────────────────────────────────────────────────
# Consolidar convergencia final (uno por morpho, el mejor rango disponible)
# ──────────────────────────────────────────────────────────────────────────
rank_order = ["Phylum", "Class", "Order", "Infraorder", "Infraclass", "Family", "Genus"]
def rank_rank(r): return rank_order.index(r) if r in rank_order else 99

final = (bridge.assign(_r=bridge.target_rank_bold.map(rank_rank))
                .sort_values(["morpho_taxon", "verified_in_bold", "_r"],
                             ascending=[True, False, True])
                .drop_duplicates("morpho_taxon", keep="first")
                .drop(columns="_r"))

# Estado final
def final_status(row):
    if row.verified_in_bold:
        return "CONVERGES"
    return "NO_CONVERGENCE_YET"
final["convergence_status"] = final.apply(final_status, axis=1)
final.to_csv(OUT / "stage3_final_convergence.tsv", sep="\t", index=False)

# ──────────────────────────────────────────────────────────────────────────
# Resumen legible
# ──────────────────────────────────────────────────────────────────────────
n_conv = (final.convergence_status == "CONVERGES").sum()
n_absent = (final.convergence_status != "CONVERGES").sum()

lines = [
    "# Stage 3 — Puente WoRMS y convergencia final",
    "",
    f"- Taxa morfológicos: **{len(CURATED)}**",
    f"- Con mapping curado que converge con BOLD local: **{n_conv}**",
    f"- Sin convergencia local (aún no representados en los 261 ASVs): **{n_absent}**",
    "",
    "## Convergencia final",
    "",
    "| Morfológico | WoRMS name | Rango objetivo | Convergencia | ASVs BOLD | Nota |",
    "|---|---|---|---|---|---|",
]
for _, r in final.iterrows():
    lines.append(
        f"| {r['morpho_taxon']} | {r['worms_name']} | {r['target_rank_bold']} | "
        f"{r['convergence_status']} | {r['n_asvs_bold']} | {r['note']} |"
    )
lines += [
    "",
    "## Comentarios interpretativos",
    "",
    "- **Convergencia esperable en Class/Order** para copépodos, cladóceros, "
    "poliquetos y cnidarios: WoRMS y BOLD comparten estos niveles Linneanos.",
    "- **Estadios larvarios (Nauplio, Zoea, Megalopa, Larva_veliger, "
    "Larva_briozoa, Huevos_pez)** no tienen registro directo en BOLD; se "
    "reconcilian a nivel de Class del adulto (Copepoda, Decapoda, "
    "Gastropoda, Bryozoa, Teleostei).",
    "- **Nombres subsumidos** (Poecilostomatoida → Cyclopoida, Mysidacea → "
    "Mysida) requieren el mapping curado; el grep léxico solo no basta.",
    "- **Grupos raros o ausentes en el marker COI Leray-XT** (Doliolida, "
    "Appendicularia, Nematoda, Euphausiidae) probablemente reflejan sesgo "
    "de amplificación / sesgo de la DB BOLD, no ausencia biológica real.",
]
(OUT / "stage3_summary.md").write_text("\n".join(lines))
print("\n".join(lines[:12]))
print()
print(final[["morpho_taxon","worms_name","target_rank_bold","convergence_status","n_asvs_bold"]].to_string(index=False))
