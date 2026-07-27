#!/usr/bin/env python3
"""
Stage 3b — Comparación cuantitativa morfología ↔ metabarcoding en las 8
estaciones EN1–EN8 usando el puente WoRMS→BOLD de Stage 3.

Pipeline:
  1. Colapsar reads BOLD a nivel Class (mejor rango de convergencia).
  2. Colapsar conteos morfológicos a los mismos targets Class (agregando
     estadios larvarios al Class del adulto).
  3. Restringir a los ASV L1..L8 (que corresponden a estaciones EN1..EN8) y
     alinear métricas relativas.
  4. Correlacionar rangos de abundancia relativa por Class y producir figuras.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

sns.set_theme(style="whitegrid", context="talk")

LANCO = Path("/sessions/gracious-sleepy-edison/mnt/LANCO")
IN    = LANCO / "results/MORPHO_METABAR"
FIG   = IN / "figures"
FIG.mkdir(exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────
# Datos morfología (long) + puente WoRMS
# ──────────────────────────────────────────────────────────────────────────
morpho = pd.read_csv(IN / "tables/morpho_long.tsv", sep="\t").rename(columns={"Taxon_morpho": "morpho_taxon"})
bridge = pd.read_csv(IN / "reports/stage3_final_convergence.tsv", sep="\t")
bridge = bridge.rename(columns={"worms_name": "Target"})[["morpho_taxon", "Target", "target_rank_bold"]]
morpho = morpho.merge(bridge, on="morpho_taxon")

# Colapsar morfología al Target (mismo nombre → suma)
morpho_class = (morpho.groupby(["Estacion", "Target"], observed=True)
                       ["absoluta"].sum().reset_index()
                       .rename(columns={"absoluta": "morpho_abs"}))

# ──────────────────────────────────────────────────────────────────────────
# Datos metabarcoding (L1..L8 = EN1..EN8, ver metadata)
# ──────────────────────────────────────────────────────────────────────────
meta = pd.read_csv(LANCO / "RUN29bis.xlsx - Metadata_COI_LANCO.tsv", sep="\t")
meta = meta.rename(columns={"ID_library":"sample_id","Estación":"Station","Year":"Year"})
meta["Station"] = meta["Station"].replace({"NE2": "EN2"})
en24 = meta[(meta.Year == 24) & (meta.Station.str.startswith("EN"))]
en24 = en24[["sample_id","Station"]]

# Class-level counts (48 ASVs Metazoa) — usamos Class que es donde más taxa convergen
cls_counts = pd.read_csv(LANCO / "results/EDA_TAX_RUN29/tables/counts_class.tsv", sep="\t").set_index("sample_id")
cls_en24 = cls_counts.loc[cls_counts.index.isin(en24.sample_id)]
cls_en24 = cls_en24.merge(en24, left_index=True, right_on="sample_id").drop(columns="sample_id").set_index("Station")
# Long
meta_class = cls_en24.stack().rename("bold_abs").reset_index()
meta_class.columns = ["Estacion", "Target", "bold_abs"]

# También intentar Order para taxa que están en Order (Calanoida, etc.)
ord_counts = pd.read_csv(LANCO / "results/EDA_TAX_RUN29/tables/counts_order.tsv", sep="\t").set_index("sample_id")
ord_en24 = ord_counts.loc[ord_counts.index.isin(en24.sample_id)]
ord_en24 = ord_en24.merge(en24, left_index=True, right_on="sample_id").drop(columns="sample_id").set_index("Station")
meta_order = ord_en24.stack().rename("bold_abs").reset_index()
meta_order.columns = ["Estacion", "Target", "bold_abs"]

meta_all = pd.concat([meta_class.assign(bold_rank="Class"),
                      meta_order.assign(bold_rank="Order")], ignore_index=True)

# ──────────────────────────────────────────────────────────────────────────
# Merge morfología ↔ metabarcoding
# ──────────────────────────────────────────────────────────────────────────
joined = morpho_class.merge(meta_all, on=["Estacion", "Target"], how="inner")
# Relativizar en cada matriz
joined["morpho_rel"] = joined.groupby("Estacion")["morpho_abs"].transform(lambda x: x / x.sum() if x.sum() else 0)
joined["bold_rel"]   = joined.groupby(["Estacion","bold_rank"])["bold_abs"].transform(lambda x: x / x.sum() if x.sum() else 0)
joined.to_csv(IN / "tables/joined_morpho_bold_EN.tsv", sep="\t", index=False)

# ──────────────────────────────────────────────────────────────────────────
# Correlación de rangos + scatter log-log
# ──────────────────────────────────────────────────────────────────────────
corr_rows = []
for rank in ["Class", "Order"]:
    sub = joined[joined.bold_rank == rank]
    if len(sub) < 3:
        continue
    rho, p = stats.spearmanr(sub["morpho_rel"], sub["bold_rel"])
    r, p_r = stats.pearsonr(np.log10(sub["morpho_rel"]+1e-4), np.log10(sub["bold_rel"]+1e-4))
    corr_rows.append({"bold_rank": rank, "n": len(sub), "spearman_rho": rho, "spearman_p": p,
                      "pearson_log_r": r, "pearson_log_p": p_r})
corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(IN / "tables/correlation_morpho_bold.tsv", sep="\t", index=False)

fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
for ax, rank in zip(axes, ["Class", "Order"]):
    sub = joined[joined.bold_rank == rank].copy()
    ax.scatter(sub["morpho_rel"]+1e-4, sub["bold_rel"]+1e-4,
               c={"EN1":"#1b9e77","EN2":"#d95f02","EN3":"#7570b3","EN4":"#e7298a",
                  "EN5":"#66a61e","EN6":"#e6ab02","EN7":"#a6761d","EN8":"#666666"}[sub["Estacion"].iloc[0]] if len(sub)<2 else "#333",
               s=80, edgecolor="black", alpha=0.75)
    for _, r in sub.iterrows():
        color = {"EN1":"#1b9e77","EN2":"#d95f02","EN3":"#7570b3","EN4":"#e7298a",
                 "EN5":"#66a61e","EN6":"#e6ab02","EN7":"#a6761d","EN8":"#666666"}[r.Estacion]
        ax.scatter(r["morpho_rel"]+1e-4, r["bold_rel"]+1e-4, c=color, s=80,
                   edgecolor="black", alpha=0.85)
        if r["morpho_rel"] > 0.05 or r["bold_rel"] > 0.05:
            ax.annotate(r["Target"], (r["morpho_rel"]+1e-4, r["bold_rel"]+1e-4),
                        fontsize=7, alpha=0.7)
    row = corr_df.query("bold_rank == @rank")
    if len(row):
        ax.set_title(f"{rank} · Spearman ρ={row['spearman_rho'].iloc[0]:.2f} "
                     f"(p={row['spearman_p'].iloc[0]:.3f}) · n={row['n'].iloc[0]}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Abundancia relativa morfología")
    ax.set_ylabel("Abundancia relativa metabarcoding")
    ax.plot([1e-4, 1], [1e-4, 1], "--", color="grey", alpha=0.5)
plt.tight_layout(); plt.savefig(FIG / "01_correlation_scatter.png", dpi=150); plt.close()

# ──────────────────────────────────────────────────────────────────────────
# Heatmap comparativo (morfología en % vs metabarcoding en %)
# ──────────────────────────────────────────────────────────────────────────
def pivot_pct(sub, val_col):
    p = sub.pivot(index="Estacion", columns="Target", values=val_col)
    return p.reindex([f"EN{i}" for i in range(1,9)])

class_sub = joined[joined.bold_rank == "Class"]
mp = pivot_pct(class_sub, "morpho_rel").fillna(0)
bp = pivot_pct(class_sub, "bold_rel").fillna(0)

fig, axes = plt.subplots(1, 2, figsize=(15, 6))
sns.heatmap(mp*100, cmap="Greens", ax=axes[0], annot=True, fmt=".0f",
            cbar_kws={"label": "% morfología"})
axes[0].set_title("Morfología (% dentro de estación)")
sns.heatmap(bp*100, cmap="Purples", ax=axes[1], annot=True, fmt=".0f",
            cbar_kws={"label": "% metabarcoding"})
axes[1].set_title("Metabarcoding COI (% dentro de estación, nivel Class)")
plt.tight_layout(); plt.savefig(FIG / "02_heatmap_side_by_side.png", dpi=150); plt.close()

# ──────────────────────────────────────────────────────────────────────────
# Barras apiladas dobles: por estación
# ──────────────────────────────────────────────────────────────────────────
targets_sorted = mp.sum(axis=0).sort_values(ascending=False).index.tolist()[:8]
palette = dict(zip(targets_sorted, sns.color_palette("tab20", len(targets_sorted))))

fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
for ax, matrix, title in [(axes[0], mp, "Morfología (formalina)"),
                          (axes[1], bp, "Metabarcoding COI-Class")]:
    bottom = np.zeros(len(matrix))
    for t in targets_sorted:
        vals = matrix[t].values if t in matrix.columns else np.zeros(len(matrix))
        ax.bar(matrix.index, vals*100, bottom=bottom, label=t, color=palette[t],
               edgecolor="white")
        bottom += vals*100
    ax.set_ylim(0, 100); ax.set_ylabel("%")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
plt.tight_layout(); plt.savefig(FIG / "03_barplot_paired.png", dpi=150); plt.close()

print("\n=== Correlación por rango ===")
print(corr_df.to_string(index=False))
print(f"\n=== Targets Class comparados: {len(mp.columns)} ===")
print(mp.columns.tolist())
