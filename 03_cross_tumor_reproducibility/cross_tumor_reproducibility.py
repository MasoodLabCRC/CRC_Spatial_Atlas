import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from itertools import combinations

# -------------------------------------------------------
# USER SETTINGS
# -------------------------------------------------------

SSGSEA_PATH  = "CRC_Tumor_Core_modules.csv"       
METADATA_PATH = "Tumor_region_meta.csv"         

#Exclude Adjacent normal samples 
SAMPLES_TO_EXCLUDE = [
    'S21_10041_B1',
    'S21_11953_A1',
]

# The region column name and expected values
REGION_COL = 'Tumor_region'
# Order regions from core to edge for display
REGION_ORDER = ['Core', 'Inner margin', 'Outer margin', 'Edge']



# -------------------------------------------------------
# LOAD & MERGE
# -------------------------------------------------------
print("Loading data...")
ssgsea = pd.read_csv(SSGSEA_PATH, sep=',')
meta   = pd.read_csv(METADATA_PATH, sep=',')

# Standardise spot ID column name — both files use barcode as first col or cell_id
# Rename to common key 'spot_id'
if 'cell_id' in meta.columns:
    meta = meta.rename(columns={'cell_id': 'barcode'})

# Merge on barcode
df = pd.merge(ssgsea, meta[['barcode', 'orig.ident', REGION_COL]],
              on='barcode', how='inner')

print(f"Merged: {len(df)} spots")
print(f"Regions found: {df[REGION_COL].unique()}")

# Exclude samples
df = df[~df['sample_id'].isin(SAMPLES_TO_EXCLUDE)]
print(f"After exclusion: {df['sample_id'].nunique()} samples, {len(df)} spots")

# Module columns
module_cols = sorted(
    [c for c in df.columns if c.startswith('Module_')],
    key=lambda x: int(x.split('_')[1])
)
print(f"Modules: {len(module_cols)}")

# -------------------------------------------------------
# STEP 1: Compute per-tumor per-region mean ssGSEA score
# -------------------------------------------------------
region_profiles = df.groupby(['sample_id', REGION_COL])[module_cols].mean().reset_index()
print(f"\nRegion profiles shape: {region_profiles.shape}")
print(region_profiles.head(3))

# -------------------------------------------------------
# STEP 2: For each module, compute pairwise Spearman
# correlation of region profiles across all tumor pairs
# -------------------------------------------------------
print("\nComputing pairwise Spearman correlations across tumor pairs...")

tumors = df['sample_id'].unique()
tumor_pairs = list(combinations(tumors, 2))
print(f"Number of tumor pairs: {len(tumor_pairs)}")

pairwise_r   = {m: [] for m in module_cols}
pairwise_p   = {m: [] for m in module_cols}

for t1, t2 in tumor_pairs:
    prof1 = region_profiles[region_profiles['sample_id'] == t1]\
            .set_index(REGION_COL)[module_cols]
    prof2 = region_profiles[region_profiles['sample_id'] == t2]\
            .set_index(REGION_COL)[module_cols]

    # Only use regions present in BOTH tumors
    common_regions = prof1.index.intersection(prof2.index)

    if len(common_regions) < 3:
        # Need at least 3 regions to compute Spearman
        continue

    for mod in module_cols:
        v1 = prof1.loc[common_regions, mod].values
        v2 = prof2.loc[common_regions, mod].values
        if np.std(v1) == 0 or np.std(v2) == 0:
            continue
        r, p = spearmanr(v1, v2)
        pairwise_r[mod].append(r)
        pairwise_p[mod].append(p)

# -------------------------------------------------------
# STEP 3: Summarise — mean pairwise r per module
# Fisher z-transform before averaging
# -------------------------------------------------------
def fisher_z(r):
    r = np.clip(r, -0.9999, 0.9999)
    return 0.5 * np.log((1 + r) / (1 - r))

def fisher_z_inv(z):
    return (np.exp(2 * z) - 1) / (np.exp(2 * z) + 1)

results = []
for mod in module_cols:
    rs = np.array(pairwise_r[mod])
    if len(rs) == 0:
        continue
    zs = fisher_z(rs)
    mean_z  = np.mean(zs)
    se_z    = np.std(zs, ddof=1) / np.sqrt(len(zs))
    mean_r  = fisher_z_inv(mean_z)
    # One-sample t-test: is mean Z > 0?
    from scipy.stats import ttest_1samp
    t, p = ttest_1samp(zs, popmean=0)
    results.append({
        'module':    mod,
        'mean_r':    mean_r,
        'mean_z':    mean_z,
        'se_z':      se_z,
        'n_pairs':   len(rs),
        'p_value':   p
    })

results_df = pd.DataFrame(results)

# BH correction
_, fdr, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
results_df['FDR'] = fdr
results_df = results_df.sort_values('mean_r', ascending=False)

print("\n--- Spearman reproducibility results ---")
print(results_df[['module', 'mean_r', 'n_pairs', 'p_value', 'FDR']].to_string(index=False))
results_df.to_csv('spearman_reproducibility.csv', index=False)

# -------------------------------------------------------
# STEP 4: PLOT 1 — Bar chart of mean pairwise Spearman r
# per module, colored by category
# -------------------------------------------------------
module_categories = {
    'Module_11': 'Canonical',   'Module_21': 'Canonical',
    'Module_1':  'Canonical',   'Module_26': 'Canonical',
    'Module_12': 'Canonical',   'Module_16': 'Canonical',
    'Module_27': 'Canonical',
    'Module_2':  'Regenerative/Fetal', 'Module_3':  'Regenerative/Fetal',
    'Module_10': 'Regenerative/Fetal', 'Module_28': 'Regenerative/Fetal',
    'Module_4':  'Non-canonical', 'Module_8':  'Non-canonical',
    'Module_22': 'Non-canonical',
    'Module_17': 'Invasive/Adaptive', 'Module_29': 'Invasive/Adaptive',
    'Module_9':  'Invasive/Adaptive', 'Module_18': 'Invasive/Adaptive',
    'Module_14': 'Metabolic', 'Module_15': 'Metabolic',
    'Module_19': 'Metabolic', 'Module_24': 'Metabolic',
    'Module_25': 'Metabolic',
    'Module_20': 'Proliferative',
    'Module_5':  'Microenvironmental', 'Module_6':  'Microenvironmental',
    'Module_7':  'Microenvironmental', 'Module_13': 'Microenvironmental',
    'Module_23': 'Microenvironmental', 'Module_30': 'Microenvironmental',
}

category_colors = {
    'Canonical':           '#2171b5',
    'Regenerative/Fetal':  '#cb181d',
    'Non-canonical':       '#6a51a3',
    'Invasive/Adaptive':   '#d94801',
    'Metabolic':           '#238b45',
    'Proliferative':       '#7a0177',
    'Microenvironmental':  '#016c59',
}

# Sort results by category then module number for grouped display
category_order = ['Canonical', 'Regenerative/Fetal', 'Non-canonical',
                  'Invasive/Adaptive', 'Metabolic', 'Proliferative',
                  'Microenvironmental']

results_df['category'] = results_df['module'].map(module_categories)
results_df['mod_num']  = results_df['module'].apply(lambda x: int(x.split('_')[1]))
results_df['cat_order'] = results_df['category'].apply(
    lambda x: category_order.index(x) if x in category_order else 99
)
results_df_plot = results_df.sort_values(['cat_order', 'mod_num'])

fig, ax = plt.subplots(figsize=(14, 5))

colors = [category_colors.get(results_df_plot.loc[i, 'category'], '#999999')
          for i in results_df_plot.index]

bars = ax.bar(
    range(len(results_df_plot)),
    results_df_plot['mean_r'],
    color=colors,
    edgecolor='white',
    linewidth=0.5,
    width=0.75
)

# Error bars (SE of Fisher z, back-transformed approximately)
ax.errorbar(
    range(len(results_df_plot)),
    results_df_plot['mean_r'],
    yerr=results_df_plot['se_z'],   # approximate SE
    fmt='none', color='black',
    capsize=2, linewidth=0.8, alpha=0.6
)

# FDR significance markers
for i, (_, row) in enumerate(results_df_plot.iterrows()):
    if row['FDR'] < 0.001:
        ax.text(i, row['mean_r'] + row['se_z'] + 0.02, '***',
                ha='center', va='bottom', fontsize=7)
    elif row['FDR'] < 0.01:
        ax.text(i, row['mean_r'] + row['se_z'] + 0.02, '**',
                ha='center', va='bottom', fontsize=7)
    elif row['FDR'] < 0.05:
        ax.text(i, row['mean_r'] + row['se_z'] + 0.02, '*',
                ha='center', va='bottom', fontsize=7)

# X axis labels
short_labels = [m.replace('Module_', 'M') for m in results_df_plot['module']]
ax.set_xticks(range(len(results_df_plot)))
ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)

# Category background shading
cat_positions = {}
for cat in category_order:
    idxs = [i for i, m in enumerate(results_df_plot['module'])
            if module_categories.get(m) == cat]
    if idxs:
        cat_positions[cat] = (min(idxs), max(idxs))

for cat, (start, end) in cat_positions.items():
    ax.axvspan(start - 0.4, end + 0.4,
               color=category_colors[cat], alpha=0.07, zorder=0)
    ax.text((start + end) / 2, -0.07,
            cat.replace('/', '/\n'), ha='center', va='top',
            fontsize=6.5, color=category_colors[cat],
            transform=ax.get_xaxis_transform())

ax.axhline(0, color='black', linewidth=0.8, linestyle='-')
ax.set_ylabel('Mean Pairwise Spearman r\n(region profiles across tumor pairs)', fontsize=10)
ax.set_title('Cross-tumor reproducibility of spatial module region profiles\n'
             '(Spearman correlation of Core/Inner/Outer/Edge mean scores)',
             fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(-0.5, len(results_df_plot) - 0.5)

plt.tight_layout()
plt.savefig('spearman_region_reproducibility.pdf',
            dpi=300, bbox_inches='tight')
plt.savefig('spearman_region_reproducibility.png',
            dpi=300, bbox_inches='tight')
print("\nSaved bar chart.")


print("\nAll done. Output files:")
print("  spearman_reproducibility.csv")
print("  spearman_region_reproducibility.pdf/.png")

