#!/usr/bin/env python3
"""
EDA — LANCO/RUN29 COI metabarcoding (38 muestras L*)
Inputs: asv_table.tsv (498 ASVs × 38 samples, pipeline PE),
        track_reads.tsv, Metadata_COI_LANCO.tsv
Outputs: figures/, tables/, eda_summary.json
"""
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import squareform
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from sklearn.manifold import MDS

warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid", context="talk")
PAL = {"EN": "#1b9e77", "H": "#d95f02", "24": "#7570b3", "25": "#e7298a"}

UPLOADS = Path("/sessions/jolly-determined-darwin/mnt/uploads")
OUT = Path("/sessions/jolly-determined-darwin/mnt/outputs/eda_lanco")
FIG = OUT / "figures"
TAB = OUT / "tables"

# ---------------------------------------------------------------------------
# 1. Carga y armonización ----------------------------------------------------
# ---------------------------------------------------------------------------
asv = pd.read_csv(UPLOADS / "asv_table.tsv", sep="\t").set_index("sample_id")
track = pd.read_csv(UPLOADS / "track_reads.tsv", sep="\t")
meta = pd.read_csv(UPLOADS / "RUN29bis.xlsx - Metadata_COI_LANCO.tsv", sep="\t")
meta = meta.rename(columns={"ID_library": "sample_id", "Estación": "Station",
                            "Year": "Year", "ID_ref": "ID_ref"})
# Corregir typo en metadata: NE2 -> EN2 (estación EN identificada por ID_ref)
meta["Station"] = meta["Station"].replace({"NE2": "EN2"})
meta["StationType"] = meta["Station"].str.extract(r"([A-Z]+)")
meta["Year"] = meta["Year"].astype(str)
meta["YearStation"] = meta["Year"] + "_" + meta["StationType"]

meta = meta[["sample_id", "ID_ref", "Station", "StationType", "Year", "YearStation"]]
asv = asv.loc[meta["sample_id"]]                    # asegurar mismo orden

summary = {"n_samples": len(meta), "n_asvs_raw": asv.shape[1]}

# ---------------------------------------------------------------------------
# 2. Tracking de reads -------------------------------------------------------
# ---------------------------------------------------------------------------
track_long = track.melt(id_vars="sample_id",
                        value_vars=["input", "q20_filtered", "denoised_F",
                                    "denoised_R", "merged", "nonchim_filt"],
                        var_name="step", value_name="reads")
step_order = ["input", "q20_filtered", "denoised_F", "denoised_R",
              "merged", "nonchim_filt"]
track_long["step"] = pd.Categorical(track_long["step"], step_order)

fig, ax = plt.subplots(1, 2, figsize=(15, 6))
sns.boxplot(data=track_long, x="step", y="reads", ax=ax[0], color="#88a0a8")
ax[0].set_yscale("log"); ax[0].set_title("Reads por etapa (log)")
ax[0].tick_params(axis="x", rotation=30)

sns.barplot(data=track.sort_values("nonchim_filt"), x="sample_id",
            y="nonchim_filt", ax=ax[1], color="#1f77b4")
ax[1].set_title("Reads no-quiméricos por muestra")
ax[1].tick_params(axis="x", rotation=90); ax[1].set_yscale("log")
plt.tight_layout(); plt.savefig(FIG / "01_read_tracking.png", dpi=150)
plt.close()

# % retención
track["pct_retained"] = 100 * track["nonchim_filt"] / track["input"]
summary["retention_median_pct"] = float(track["pct_retained"].median())
summary["retention_min_pct"] = float(track["pct_retained"].min())
summary["retention_max_pct"] = float(track["pct_retained"].max())
summary["L33_reads_input"] = int(track.loc[track.sample_id == "L33", "input"].iloc[0])
summary["L33_reads_final"] = int(track.loc[track.sample_id == "L33", "nonchim_filt"].iloc[0])
summary["median_reads_final"] = int(track["nonchim_filt"].median())

# ---------------------------------------------------------------------------
# 3. Métricas globales ASV ---------------------------------------------------
# ---------------------------------------------------------------------------
asv_counts = asv.copy()
asv_prevalence = (asv_counts > 0).sum(axis=0)
asv_total_abund = asv_counts.sum(axis=0)
asv_stats = pd.DataFrame({"prevalence": asv_prevalence,
                          "total_reads": asv_total_abund})
asv_stats["singleton_sample"] = asv_stats["prevalence"] == 1
summary["n_asvs_single_sample"] = int(asv_stats["singleton_sample"].sum())
summary["n_asvs_core_50pct"] = int((asv_prevalence >= 19).sum())
asv_stats.to_csv(TAB / "asv_global_stats.tsv", sep="\t")

# Rank-abundance + prevalencia
fig, ax = plt.subplots(1, 2, figsize=(15, 5))
ax[0].plot(np.arange(1, len(asv_total_abund) + 1),
           asv_total_abund.sort_values(ascending=False).values, color="#333")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_title("Curva rank-abundance (global)")
ax[0].set_xlabel("Rango ASV"); ax[0].set_ylabel("Reads totales")

sns.histplot(asv_prevalence, bins=range(1, 40), ax=ax[1], color="#1b9e77")
ax[1].set_title("Distribución de prevalencia (muestras donde aparece cada ASV)")
ax[1].set_xlabel("Nº de muestras"); ax[1].set_ylabel("Nº de ASVs")
plt.tight_layout(); plt.savefig(FIG / "02_asv_global_stats.png", dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 4. Normalización: rarefacción al mínimo ≥ 5 000 ----------------------------
# ---------------------------------------------------------------------------
sample_sums = asv.sum(axis=1)
rar_depth = int(max(5000, sample_sums.min()))      # asegurar piso de 5K
summary["rarefaction_depth"] = rar_depth

rng = np.random.default_rng(42)
def rarefy(row, depth):
    if row.sum() < depth:
        return row.copy()
    p = row.values / row.values.sum()
    drawn = rng.multinomial(depth, p)
    return pd.Series(drawn, index=row.index)

asv_rar = asv.apply(lambda r: rarefy(r, rar_depth), axis=1)
asv_rar.to_csv(TAB / "asv_rarefied.tsv", sep="\t")

# Proporcional (alternativa para beta)
asv_rel = asv.div(asv.sum(axis=1), axis=0)
asv_rel.to_csv(TAB / "asv_relative.tsv", sep="\t")

# ---------------------------------------------------------------------------
# 5. Diversidad ALFA ---------------------------------------------------------
# ---------------------------------------------------------------------------
counts = asv_rar.values.astype(int)
ids = asv_rar.index.tolist()
alpha_df = pd.DataFrame({
    "observed":   alpha_diversity("observed_otus", counts, ids=ids),
    "shannon":    alpha_diversity("shannon",       counts, ids=ids),
    "simpson":    alpha_diversity("simpson",       counts, ids=ids),
    "chao1":      alpha_diversity("chao1",         counts, ids=ids),
    "pielou_e":   alpha_diversity("pielou_e",      counts, ids=ids),
})
alpha_df.index.name = "sample_id"
alpha_df = alpha_df.merge(meta, on="sample_id", how="left")
alpha_df.to_csv(TAB / "alpha_diversity.tsv", sep="\t", index=False)

# Kruskal-Wallis por Year y StationType
def kw(df, metric, group):
    groups = [g[metric].values for _, g in df.groupby(group) if len(g) > 1]
    if len(groups) < 2:
        return np.nan, np.nan
    h, p = stats.kruskal(*groups)
    return h, p

alpha_stats = []
for m in ["observed", "shannon", "simpson", "chao1", "pielou_e"]:
    for g in ["Year", "StationType", "YearStation"]:
        h, p = kw(alpha_df, m, g)
        alpha_stats.append({"metric": m, "group": g, "H": h, "p": p})
alpha_stats_df = pd.DataFrame(alpha_stats)
alpha_stats_df.to_csv(TAB / "alpha_kruskal.tsv", sep="\t", index=False)

# Figura alfa por grupos
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
for ax, m in zip(axes.flat, ["observed", "shannon", "chao1", "pielou_e"]):
    sns.boxplot(data=alpha_df, x="YearStation", y=m, ax=ax,
                palette=["#1b9e77", "#d95f02", "#a6cee3", "#fb9a99"])
    sns.stripplot(data=alpha_df, x="YearStation", y=m, ax=ax, color="black",
                  size=4, alpha=0.6)
    ax.set_title(f"{m.title()}")
    ax.set_xlabel("")
plt.suptitle("Diversidad alfa por Year × StationType", fontsize=16)
plt.tight_layout(); plt.savefig(FIG / "03_alpha_diversity.png", dpi=150)
plt.close()

summary["alpha_observed_median"] = float(alpha_df["observed"].median())
summary["alpha_shannon_median"]  = float(alpha_df["shannon"].median())

# ---------------------------------------------------------------------------
# 6. Diversidad BETA + ordenaciones -----------------------------------------
# ---------------------------------------------------------------------------
bc = beta_diversity("braycurtis", counts, ids=ids)
ja = beta_diversity("jaccard",    (counts > 0).astype(int), ids=ids)

bc.to_data_frame().to_csv(TAB / "beta_braycurtis.tsv", sep="\t")
ja.to_data_frame().to_csv(TAB / "beta_jaccard.tsv",    sep="\t")

# PCoA Bray-Curtis
ord_bc = pcoa(bc)
exp_bc = ord_bc.proportion_explained.iloc[:3].values
coord_bc = ord_bc.samples.iloc[:, :3].copy()
coord_bc.index = ids
coord_bc = coord_bc.merge(meta, left_index=True, right_on="sample_id")

# NMDS Bray-Curtis (sklearn MDS sobre distancia, 2D)
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42,
          n_init=10, max_iter=500)
nmds_coord = mds.fit_transform(squareform(bc.condensed_form()))
coord_nmds = pd.DataFrame(nmds_coord, columns=["NMDS1", "NMDS2"], index=ids)
coord_nmds = coord_nmds.merge(meta, left_index=True, right_on="sample_id")
summary["nmds_stress"] = float(mds.stress_)

# Figuras ordenación
fig, ax = plt.subplots(1, 2, figsize=(15, 6.5))
for grp, color in zip(["EN", "H"], ["#1b9e77", "#d95f02"]):
    sub = coord_bc[coord_bc.StationType == grp]
    ax[0].scatter(sub["PC1"], sub["PC2"], s=110, c=color, edgecolor="black",
                  label=f"Tipo {grp}", alpha=0.85)
for yr, marker in zip(["24", "25"], ["o", "^"]):
    sub = coord_bc[coord_bc.Year == yr]
    ax[0].scatter(sub["PC1"], sub["PC2"], s=40, marker=marker,
                  facecolors="none", edgecolor="black", label=f"Year 20{yr}")
ax[0].set_xlabel(f"PCo1 ({exp_bc[0]*100:.1f}%)")
ax[0].set_ylabel(f"PCo2 ({exp_bc[1]*100:.1f}%)")
ax[0].set_title("PCoA Bray-Curtis")
ax[0].legend(loc="best", fontsize=10)
for _, r in coord_bc.iterrows():
    ax[0].annotate(r.sample_id, (r.PC1, r.PC2), fontsize=7, alpha=0.6)

for grp, color in zip(["EN", "H"], ["#1b9e77", "#d95f02"]):
    sub = coord_nmds[coord_nmds.StationType == grp]
    ax[1].scatter(sub["NMDS1"], sub["NMDS2"], s=110, c=color,
                  edgecolor="black", label=f"Tipo {grp}", alpha=0.85)
ax[1].set_title(f"NMDS Bray-Curtis (stress={mds.stress_:.3f})")
ax[1].set_xlabel("NMDS1"); ax[1].set_ylabel("NMDS2")
ax[1].legend(loc="best", fontsize=10)
for _, r in coord_nmds.iterrows():
    ax[1].annotate(r.sample_id, (r.NMDS1, r.NMDS2), fontsize=7, alpha=0.6)
plt.tight_layout(); plt.savefig(FIG / "04_ordination.png", dpi=150)
plt.close()

# PERMANOVA
perm_results = []
for grp in ["Year", "StationType", "YearStation"]:
    grouping = meta.set_index("sample_id").loc[ids, grp].values
    for dist, name in [(bc, "braycurtis"), (ja, "jaccard")]:
        res = permanova(dist, grouping, permutations=999)
        perm_results.append({
            "distance": name, "grouping": grp,
            "F": float(res["test statistic"]),
            "p": float(res["p-value"]),
            "n_groups": int(res["number of groups"]),
        })
perm_df = pd.DataFrame(perm_results)
perm_df.to_csv(TAB / "permanova.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# 7. Heatmap top-30 ASVs -----------------------------------------------------
# ---------------------------------------------------------------------------
top30 = asv_rel.sum(axis=0).nlargest(30).index
hm = asv_rel[top30].T
sample_order = coord_bc.sort_values(["StationType", "Year"]).sample_id.tolist()
hm = hm[sample_order]
fig, ax = plt.subplots(figsize=(14, 9))
sns.heatmap(np.log10(hm + 1e-4), cmap="rocket_r", ax=ax,
            cbar_kws={"label": "log10(abundancia relativa)"})
ax.set_title("Top-30 ASVs (abundancia relativa rarefied)")
plt.tight_layout(); plt.savefig(FIG / "05_heatmap_top30.png", dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 8. Persistir resumen JSON --------------------------------------------------
# ---------------------------------------------------------------------------
summary["permanova"] = perm_results
summary["alpha_kruskal"] = alpha_stats
with open(OUT / "eda_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=str)

print(json.dumps(summary, indent=2, default=str))
