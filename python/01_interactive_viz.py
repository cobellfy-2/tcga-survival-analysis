"""
01_interactive_viz.py
Loads preprocessed results from R and produces interactive Plotly figures.
Outputs: results/figures/km_interactive.html, volcano_interactive.html
"""

import os
import sys
import subprocess

# Auto-install missing packages
def require(pkg, import_name=None):
    import_name = import_name or pkg
    try:
        __import__(import_name)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])

for p in [("pyreadr", "pyreadr"), ("plotly", "plotly"), ("pandas", "pandas"),
          ("numpy", "numpy"), ("lifelines", "lifelines"), ("kaleido", "kaleido")]:
    require(*p)

import pyreadr
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from lifelines import KaplanMeierFitter
from pathlib import Path

RESULTS_DIR = Path("results")
FIGURES_DIR = Path("results/figures")
PROCESSED_DIR = Path("data/processed")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# --- Load data (written by R with saveRDS, read via pyreadr) ---
# pyreadr reads .rds files
survival_rds  = pyreadr.read_r(str(PROCESSED_DIR / "survival_df.rds"))[None]
gene_cox_rds  = pyreadr.read_r(str(RESULTS_DIR / "gene_cox_results.rds"))[None]

# ── Interactive Volcano ──────────────────────────────────────────────────────
df_cox = gene_cox_rds.copy()
df_cox["neg_log10p"] = -np.log10(df_cox["padj"].clip(lower=1e-300))
df_cox["significance"] = "Not significant"
df_cox.loc[(df_cox["padj"] < 0.05) & (df_cox["log2HR"] > 0), "significance"] = "High = worse OS"
df_cox.loc[(df_cox["padj"] < 0.05) & (df_cox["log2HR"] < 0), "significance"] = "High = better OS"

color_map = {
    "High = worse OS":  "#F44336",
    "High = better OS": "#2196F3",
    "Not significant":  "#BDBDBD",
}

fig_volcano = px.scatter(
    df_cox,
    x="log2HR",
    y="neg_log10p",
    color="significance",
    color_discrete_map=color_map,
    hover_name="gene",
    hover_data={"hr": ":.3f", "padj": ":.2e", "log2HR": ":.3f"},
    labels={"log2HR": "log2 Hazard Ratio", "neg_log10p": "-log10(adj. p-value)"},
    title="TCGA-BRCA: Gene-level Survival Association (interactive)",
    opacity=0.6,
)
fig_volcano.add_vline(x=0, line_dash="dash", line_color="grey")
fig_volcano.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey")
fig_volcano.update_layout(legend_title_text="")

volcano_path = str(FIGURES_DIR / "volcano_interactive.html")
fig_volcano.write_html(volcano_path)
print(f"Saved: {volcano_path}")

# ── Interactive Kaplan-Meier (lifelines + plotly) ────────────────────────────
df_surv = survival_rds.copy()
df_surv["OS.time_months"] = pd.to_numeric(df_surv["OS.time"], errors="coerce") / 30.44
df_surv["OS"] = pd.to_numeric(df_surv["OS"], errors="coerce")
import re
def classify_stage(s):
    if not isinstance(s, str):
        return None
    if re.search(r"Stage IV", s):   return "IV"
    if re.search(r"Stage III", s):  return "III"
    if re.search(r"Stage II", s):   return "II"
    if re.search(r"Stage I", s):    return "I"
    return None

df_surv["stage_clean"] = df_surv["ajcc_pathologic_tumor_stage"].apply(classify_stage)
df_surv = df_surv.dropna(subset=["OS.time_months", "OS", "stage_clean"])
df_surv = df_surv[df_surv["OS.time_months"] > 0]

fig_km = go.Figure()
colors = {"I": "#2196F3", "II": "#4CAF50", "III": "#FF9800", "IV": "#F44336"}

for stage, grp in df_surv.groupby("stage_clean"):
    kmf = KaplanMeierFitter()
    kmf.fit(grp["OS.time_months"], event_observed=grp["OS"], label=stage)
    timeline = kmf.survival_function_.index
    sf       = kmf.survival_function_[stage].values
    ci_low   = kmf.confidence_interval_[f"{stage}_lower_0.95"].values
    ci_high  = kmf.confidence_interval_[f"{stage}_upper_0.95"].values

    fig_km.add_trace(go.Scatter(
        x=np.concatenate([timeline, timeline[::-1]]),
        y=np.concatenate([ci_high, ci_low[::-1]]),
        fill="toself", fillcolor=colors[stage],
        opacity=0.15, line=dict(width=0),
        showlegend=False, name=f"{stage} CI",
    ))
    fig_km.add_trace(go.Scatter(
        x=timeline, y=sf,
        mode="lines", name=f"Stage {stage}",
        line=dict(color=colors[stage], width=2),
    ))

fig_km.update_layout(
    title="TCGA-BRCA: Overall Survival by Tumor Stage (interactive)",
    xaxis_title="Time (months)",
    yaxis_title="Survival Probability",
    yaxis=dict(range=[0, 1.05]),
    legend_title="Stage",
    template="plotly_white",
)

km_path = str(FIGURES_DIR / "km_interactive.html")
fig_km.write_html(km_path)
print(f"Saved: {km_path}")

print("\nInteractive figures complete. Next: report/report.qmd")
