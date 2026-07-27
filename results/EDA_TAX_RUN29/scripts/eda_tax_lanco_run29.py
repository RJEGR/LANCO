#!/usr/bin/env python3
"""
EDA taxonómico — LANCO/RUN29 COI SE pipeline (261 ASVs × 38 muestras L)
Inputs:
    asv_table.tsv              (261 ASVs × 38 muestras, SE)
    tax_BOLD.tsv               (Feature ID, Taxon k__;p__;c__;o__;f__;g__;s__;t__, Consensus)
    track_reads.tsv            (SE tracking)
    RUN29bis.xlsx — Metadata_COI_LANCO.tsv
Outputs: figures/, tables/, eda_tax_summary.json
"""
import json
import re
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
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid", context="talk")

UPLOADS = Path("/sessions/trusting-lucid-cori/mnt/uploads")
META_TSV = Path("/sessions/trusting-lucid-cori/mnt/LANCO/RUN29bis.xlsx - Metadata_COI_LANCO.tsv")
OUT = Path("/sessions/trusting-lucid-cori/mnt/outputs/eda_tax")
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

RANK_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "BOLD_BIN"]
RANK_PREFIX = ["k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__"]

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 · Cargar y parsear tax_BOLD ---------------------------------------
# ═══════════════════════════════════════════════════════════════════════════
tax_raw = pd.read_csv(UPLOADS / "tax_BOLD.tsv", sep="\t")
tax_raw.columns = ["ASV_ID", "Taxon", "Consensus"]

def parse_taxon(s):
    """Devuelve dict con 8 rangos; NaN si no está presente"""
    out = {r: np.nan for r in RANK_LEVELS}
    if pd.isna(s) or s.strip() in ("Unassigned", ""):
        return out
    parts = [p.strip() for p in str(s).split(";")]
    for p in parts:
        for pref, rank in zip(RANK_PREFIX, RANK_LEVELS):
            if p.startswith(pref):
                val = p[len(pref):].strip()
                out[rank] = val if val else np.nan
                break
    return out

parsed = tax_raw["Taxon"].apply(parse_taxon).apply(pd.Series)
tax = pd.concat([tax_raw[["ASV_ID", "Consensus"]], parsed], axis=1).set_index("ASV_ID")
tax["Kingdom"] = tax["Kingdom"].fillna("Unassigned")
tax.to_csv(TAB / "tax_parsed.tsv", sep="\t")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1b · Cargar counts + metadata ---------------------------------------
# ═══════════════════════════════════════════════════════════════════════════
asv = pd.read_csv(UPLOADS / "asv_table.tsv", sep="\t").set_index("sample_id")
track = pd.read_csv(UPLOADS / "track_reads.tsv", sep="\t")
meta = pd.read_csv(META_TSV, sep="\t")
meta = meta.rename(columns={"ID_library": "sample_id", "Estación": "Station", "Year": "Year"})
meta["Station"] = meta["Station"].replace({"NE2": "EN2"})
meta["StationType"] = meta["Station"].str.extract(r"([A-Z]+)")
meta["Year"] = meta["Year"].astype(str)
meta["YearStation"] = meta["Year"] + "_" + meta["StationType"]
asv = asv.loc[meta["sample_id"]]

summary = {
    "n_samples": int(len(meta)),
    "n_asvs_raw": int(asv.shape[1]),
    "kingdom_counts": tax["Kingdom"].value_counts().to_dict(),
    "consensus_median": float(tax["Consensus"].median()),
    "consensus_min": float(tax["Consensus"].min()),
    "phyla_animalia": tax.loc[tax.Kingdom == "Animalia", "Phylum"].value_counts().to_dict(),
}

# ═══════════════════════════════════════════════════════════════════════════
# STEP 2 · Filtro Metazoa + consensus + prevalencia ------------------------
# ═══════════════════════════════════════════════════════════════════════════
CONS_MIN = 0.7
PREVALENCE_MIN = 2

metazoa_ids = tax[(tax["Kingdom"] == "Animalia") &
                  (tax["Consensus"] >= CONS_MIN)].index
asv_meta = asv.loc[:, asv.columns.intersection(metazoa_ids)]
prev = (asv_meta > 0).sum(axis=0)
asv_meta = asv_meta.loc[:, prev >= PREVALENCE_MIN]

summary["n_asvs_metazoa_raw"] = int(len(metazoa_ids))
summary["n_asvs_metazoa_filtered"] = int(asv_meta.shape[1])
summary["reads_pre_filter"] = int(asv.values.sum())
summary["reads_post_filter"] = int(asv_meta.values.sum())
summary["pct_reads_metazoa"] = round(100 * summary["reads_post_filter"] / summary["reads_pre_filter"], 2)

tax_meta = tax.loc[asv_meta.columns].copy()
tax_meta.to_csv(TAB / "tax_metazoa_filtered.tsv", sep="\t")
asv_meta.to_csv(TAB / "asv_metazoa_counts.tsv", sep="\t")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 · Colapso a rangos taxonómicos ------------------------------------
# ═══════════════════════════════════════════════════════════════════════════
def collapse(counts_df, tax_df, rank):
    """Colapsa columnas ASV por rank agrupando y sumando; NA → 'Unclassified_{rank}'"""
    labels = tax_df[rank].fillna(f"Unclassified_{rank}")
    df = counts_df.T.copy()
    df["_lbl"] = labels.values
    grp = df.groupby("_lbl").sum().T
    grp.index.name = "sample_id"
    return grp

collapsed = {rank: collapse(asv_meta, tax_meta, rank) for rank in
             ["Phylum", "Class", "Order", "Family", "Genus"]}
for r, tbl in collapsed.items():
    tbl.to_csv(TAB / f"counts_{r.lower()}.tsv", sep="\t")
    summary[f"n_{r.lower()}"] = int(tbl.shape[1])

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 · Barplots composición Year × StationType --------------------------
# ═══════════════════════════════════════════════════════════════════════════
def stacked_barplot(counts, tax_rank, top_n, out_png, group_col="YearStation"):
    rel = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    top = rel.sum(axis=0).nlargest(top_n).index
    rel_top = rel[top].copy()
    rel_top["Other"] = 1 - rel_top.sum(axis=1)
    rel_top = rel_top.merge(meta[["sample_id", group_col, "StationType", "Year"]],
                            left_index=True, right_on="sample_id")
    rel_top = rel_top.sort_values(["StationType", "Year", "sample_id"])
    x = rel_top["sample_id"].values
    bottom = np.zeros(len(rel_top))
    palette = sns.color_palette("tab20", top_n) + [(0.6, 0.6, 0.6)]
    fig, ax = plt.subplots(figsize=(15, 7))
    for taxon, color in zip(list(top) + ["Other"], palette):
        vals = rel_top[taxon].values
        ax.bar(x, vals, bottom=bottom, label=taxon, color=color, edgecolor="white", linewidth=0.3)
        bottom += vals
    ax.set_ylabel("Abundancia relativa"); ax.set_ylim(0, 1)
    ax.set_title(f"Composición por muestra · {tax_rank} (top {top_n})")
    ax.tick_params(axis="x", rotation=90, labelsize=9)
    # Separadores entre StationType
    prev_type = None
    for i, (_, r) in enumerate(rel_top.iterrows()):
        if prev_type is not None and r["StationType"] != prev_type:
            ax.axvline(i - 0.5, color="black", linewidth=1.2, alpha=0.5)
        prev_type = r["StationType"]
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, ncol=1)
    plt.tight_layout(); plt.savefig(out_png, dpi=150); plt.close()

stacked_barplot(collapsed["Phylum"], "Phylum", top_n=8,  out_png=FIG / "01_barplot_phylum.png")
stacked_barplot(collapsed["Class"],  "Class",  top_n=10, out_png=FIG / "02_barplot_class.png")
stacked_barplot(collapsed["Order"],  "Order",  top_n=12, out_png=FIG / "03_barplot_order.png")
stacked_barplot(collapsed["Family"], "Family", top_n=15, out_png=FIG / "04_barplot_family.png")
stacked_barplot(collapsed["Genus"],  "Genus",  top_n=15, out_png=FIG / "05_barplot_genus.png")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 · Alfa y beta por rango taxonómico --------------------------------
# ═══════════════════════════════════════════════════════════════════════════
rar_depth = int(max(1000, asv_meta.sum(axis=1).min()))
summary["rarefaction_depth_meta"] = rar_depth
rng = np.random.default_rng(42)

def rarefy(row, depth):
    total = row.sum()
    if total < depth:
        return row.copy()
    p = row.values / total
    return pd.Series(rng.multinomial(depth, p), index=row.index)

def compute_alpha_beta(counts, level):
    rar = counts.apply(lambda r: rarefy(r, rar_depth), axis=1)
    ids = rar.index.tolist()
    arr = rar.values.astype(int)
    ad = pd.DataFrame({
        "observed": alpha_diversity("sobs",    arr, ids=ids),
        "shannon":  alpha_diversity("shannon", arr, ids=ids),
        "simpson":  alpha_diversity("simpson", arr, ids=ids),
        "chao1":    alpha_diversity("chao1",   arr, ids=ids),
        "pielou_e": alpha_diversity("pielou_e",arr, ids=ids),
    })
    ad["level"] = level
    ad = ad.merge(meta, left_index=True, right_on="sample_id")
    bc = beta_diversity("braycurtis", arr, ids=ids)
    # PERMANOVA
    perms = []
    for grp in ["Year", "StationType", "YearStation"]:
        grouping = meta.set_index("sample_id").loc[ids, grp].values
        res = permanova(bc, grouping, permutations=999)
        perms.append({"level": level, "grouping": grp,
                      "F": float(res["test statistic"]),
                      "p": float(res["p-value"])})
    return ad, bc, perms, rar

alpha_all, perm_all, dists = [], [], {}
for level, tbl in [("ASV_meta", asv_meta)] + list(collapsed.items()):
    ad, bc, perms, rar = compute_alpha_beta(tbl, level)
    alpha_all.append(ad)
    perm_all.extend(perms)
    dists[level] = bc

alpha_df = pd.concat(alpha_all, ignore_index=True)
alpha_df.to_csv(TAB / "alpha_by_level.tsv", sep="\t", index=False)
perm_df = pd.DataFrame(perm_all)
perm_df.to_csv(TAB / "permanova_by_level.tsv", sep="\t", index=False)

# Kruskal por level × métrica × grupo
kw_rows = []
for level in alpha_df["level"].unique():
    sub = alpha_df[alpha_df.level == level]
    for m in ["observed", "shannon", "simpson", "chao1", "pielou_e"]:
        for grp in ["Year", "StationType", "YearStation"]:
            groups = [g[m].dropna().values for _, g in sub.groupby(grp) if len(g) > 1]
            if len(groups) < 2:
                continue
            try:
                H, p = stats.kruskal(*groups)
            except ValueError:
                H, p = np.nan, np.nan
            kw_rows.append({"level": level, "metric": m, "group": grp, "H": H, "p": p})
kw_df = pd.DataFrame(kw_rows)
kw_df.to_csv(TAB / "alpha_kruskal_by_level.tsv", sep="\t", index=False)

# Figura alfa por rango
levels_plot = ["ASV_meta", "Genus", "Family", "Order", "Class", "Phylum"]
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=False)
for ax, lvl in zip(axes.flat, levels_plot):
    sub = alpha_df[alpha_df.level == lvl].copy()
    sub["Year"] = sub["Year"].astype(str)
    # boxplot manual (seaborn+hue bug con dodge en boxplot vacio para algún grupo)
    positions = {"EN": 0, "H": 1}
    width = 0.35
    for yr, offset, color in [("24", -width/2, "#7570b3"), ("25", +width/2, "#e7298a")]:
        data = [sub[(sub.StationType == st) & (sub.Year == yr)]["shannon"].values
                for st in ["EN", "H"]]
        xs = [positions[st] + offset for st in ["EN", "H"]]
        bp = ax.boxplot(data, positions=xs, widths=width*0.9, patch_artist=True,
                        showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(color); patch.set_alpha(0.7)
        for i, st in enumerate(["EN", "H"]):
            y = sub[(sub.StationType == st) & (sub.Year == yr)]["shannon"].values
            ax.scatter(np.full(len(y), xs[i]) + np.random.uniform(-0.05, 0.05, len(y)),
                       y, color="black", s=15, alpha=0.7, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["EN", "H"])
    ax.set_xlabel("StationType"); ax.set_ylabel("Shannon")
    ax.set_title(f"Shannon @ {lvl}")
# leyenda global
handles = [plt.Rectangle((0,0),1,1, facecolor="#7570b3", alpha=0.7, label="Year 2024"),
           plt.Rectangle((0,0),1,1, facecolor="#e7298a", alpha=0.7, label="Year 2025")]
fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.98))
plt.suptitle("Diversidad alfa (Shannon) por rango taxonómico", fontsize=16)
plt.tight_layout(); plt.savefig(FIG / "06_alpha_shannon_by_level.png", dpi=150); plt.close()

# Figura PCoA multi-level
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
for ax, lvl in zip(axes.flat, levels_plot):
    bc = dists[lvl]
    ord_ = pcoa(bc)
    exp = ord_.proportion_explained.iloc[:2].values * 100
    coord = ord_.samples.iloc[:, :2].copy()
    coord.index = bc.ids
    coord = coord.merge(meta, left_index=True, right_on="sample_id")
    for grp, marker in zip(["EN", "H"], ["o", "s"]):
        for yr, color in zip(["24", "25"], ["#1b9e77", "#d95f02"]):
            sub = coord[(coord.StationType == grp) & (coord.Year == yr)]
            ax.scatter(sub["PC1"], sub["PC2"], s=110, marker=marker, c=color,
                       edgecolor="black", label=f"{grp}·20{yr}", alpha=0.85)
    p = float(perm_df.query("level == @lvl and grouping == 'StationType'")["p"].iloc[0])
    F = float(perm_df.query("level == @lvl and grouping == 'StationType'")["F"].iloc[0])
    ax.set_title(f"{lvl}\nPERMANOVA StationType: F={F:.2f} · p={p:.3f}")
    ax.set_xlabel(f"PCo1 ({exp[0]:.1f}%)"); ax.set_ylabel(f"PCo2 ({exp[1]:.1f}%)")
    if ax == axes[0, 0]:
        ax.legend(loc="best", fontsize=8)
plt.suptitle("PCoA Bray-Curtis por rango taxonómico", fontsize=16)
plt.tight_layout(); plt.savefig(FIG / "07_pcoa_by_level.png", dpi=150); plt.close()

# ═══════════════════════════════════════════════════════════════════════════
# STEP 6 · Composición diferencial EN vs H (Wilcoxon + FDR sobre CLR) ------
# ═══════════════════════════════════════════════════════════════════════════
def clr(x, eps=0.5):
    x = x + eps
    logx = np.log(x)
    return logx.sub(logx.mean(axis=1), axis=0)

def differential(counts, level):
    """Wilcoxon 2-sample EN vs H sobre CLR, con BH-FDR"""
    rel = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    clr_df = clr(counts)
    labels = meta.set_index("sample_id").loc[counts.index, "StationType"].values
    rows = []
    for taxon in counts.columns:
        en = clr_df.loc[labels == "EN", taxon]
        h  = clr_df.loc[labels == "H",  taxon]
        if (en > 0).sum() + (h > 0).sum() < 3:  # taxa poco prevalentes
            continue
        try:
            u, p = stats.mannwhitneyu(en, h, alternative="two-sided")
        except ValueError:
            continue
        rows.append({
            "level": level, "taxon": taxon,
            "mean_rel_EN": rel.loc[labels == "EN", taxon].mean(),
            "mean_rel_H":  rel.loc[labels == "H",  taxon].mean(),
            "prev_EN": (counts.loc[labels == "EN", taxon] > 0).sum(),
            "prev_H":  (counts.loc[labels == "H",  taxon] > 0).sum(),
            "clr_EN_median": float(en.median()),
            "clr_H_median":  float(h.median()),
            "log2FC_clr": float(en.median() - h.median()),
            "U": float(u), "p": float(p),
        })
    df = pd.DataFrame(rows)
    if len(df):
        df["p_adj"] = multipletests(df["p"], method="fdr_bh")[1]
        df["signif"] = df["p_adj"] < 0.05
    return df

diff_all = pd.concat([differential(collapsed[r], r) for r in
                      ["Phylum", "Class", "Order", "Family", "Genus"]],
                     ignore_index=True)
if len(diff_all) == 0:
    diff_all = pd.DataFrame(columns=["level", "taxon", "mean_rel_EN", "mean_rel_H",
                                     "prev_EN", "prev_H", "clr_EN_median",
                                     "clr_H_median", "log2FC_clr", "U", "p",
                                     "p_adj", "signif"])
else:
    diff_all = diff_all.sort_values(["level", "p_adj", "p"])
diff_all.to_csv(TAB / "differential_EN_vs_H.tsv", sep="\t", index=False)
summary["n_signif_genera"] = int(((diff_all.level == "Genus") & diff_all.signif).sum())
summary["n_signif_family"] = int(((diff_all.level == "Family") & diff_all.signif).sum())

# Heatmap biomarcadores género (top signif por p_adj)
sig_gen = diff_all[(diff_all.level == "Genus") & diff_all.signif].head(30)
if len(sig_gen):
    rel_gen = collapsed["Genus"].div(collapsed["Genus"].sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    hm = rel_gen[sig_gen["taxon"].values].T
    sample_order = meta.sort_values(["StationType", "Year"]).sample_id.tolist()
    hm = hm[sample_order]
    fig, ax = plt.subplots(figsize=(15, max(4, 0.35 * len(hm))))
    sns.heatmap(np.log10(hm + 1e-4), cmap="rocket_r", ax=ax,
                cbar_kws={"label": "log10 abundancia relativa"})
    ax.set_title(f"Biomarcadores por género · Wilcoxon FDR<0.05 (n={len(sig_gen)})")
    ax.set_xlabel(""); ax.set_ylabel("Género")
    plt.tight_layout(); plt.savefig(FIG / "08_heatmap_biomarkers_genus.png", dpi=150); plt.close()

# Volcano por género
gen = diff_all[diff_all.level == "Genus"].copy()
if len(gen):
    fig, ax = plt.subplots(figsize=(9, 7))
    color = np.where(gen.signif & (gen.log2FC_clr > 0), "#1b9e77",
             np.where(gen.signif & (gen.log2FC_clr < 0), "#d95f02", "#bbbbbb"))
    ax.scatter(gen["log2FC_clr"], -np.log10(gen["p_adj"]), c=color, s=60,
               edgecolor="black", linewidth=0.4, alpha=0.85)
    ax.axhline(-np.log10(0.05), ls="--", c="grey")
    ax.axvline(0, ls="--", c="grey")
    for _, r in gen[gen.signif].head(12).iterrows():
        ax.annotate(r["taxon"][:20], (r["log2FC_clr"], -np.log10(r["p_adj"])),
                    fontsize=8, alpha=0.9)
    ax.set_xlabel("Δ CLR (EN − H)"); ax.set_ylabel("−log10 FDR")
    ax.set_title("Volcano · composición diferencial por género (EN vs H)")
    plt.tight_layout(); plt.savefig(FIG / "09_volcano_genus.png", dpi=150); plt.close()

# ═══════════════════════════════════════════════════════════════════════════
# Guardar summary -----------------------------------------------------------
# ═══════════════════════════════════════════════════════════════════════════
summary["permanova_summary"] = perm_all
with open(OUT / "eda_tax_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(json.dumps({k: v for k, v in summary.items() if k != "permanova_summary"}, indent=2, default=str))
print("\n--- PERMANOVA por nivel taxonómico (Bray-Curtis) ---")
print(perm_df.to_string(index=False))
print(f"\nGéneros diferencialmente abundantes EN vs H (FDR<0.05): {summary['n_signif_genera']}")
print(f"Familias diferencialmente abundantes EN vs H (FDR<0.05): {summary['n_signif_family']}")
