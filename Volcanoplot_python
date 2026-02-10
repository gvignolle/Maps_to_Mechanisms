import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import log10

# 1) Load data
fn = 'top50_features_Brain_AAN_structural_CES.xlsx'
df = pd.read_excel(fn, engine='openpyxl')

# Ensure required columns exist
required_cols = ['logFC','P.Value','adj.P.Val','Biomarkers']
assert all(c in df.columns for c in required_cols), f"Missing columns in file: {set(required_cols) - set(df.columns)}"

# 2) Compute -log10 p-value
# Guard against p-values of 0
pvals = df['P.Value'].replace(0, np.nextafter(0, 1))
df['neglog10p'] = -np.log10(pvals)

# 3) Significance & direction
alpha = 0.05

df['sig'] = df['adj.P.Val'] < alpha

def label_direction(row):
    if not row['sig']:
        return 'NS'
    return 'Up' if row['logFC'] > 0 else 'Down'

df['direction'] = df.apply(label_direction, axis=1)

# 4) Select top 5 features by smallest p-value
df_sorted = df.sort_values('P.Value', ascending=True)
labels_top = set(df_sorted.head(5)['Biomarkers'].astype(str))

# 5) Plot (1x1 aspect)
plt.figure(figsize=(6,6))  # square

# Colors (ggplot2 defaults): red/orange #F8766D and greenblue #00BFC4
colors = {
    'Up': '#F8766D',        # orangereddish
    'Down': '#00BFC4',      # greenishblue
    'NS': 'lightgray'
}

# Scatter by groups for clear legend
for grp in ['NS', 'Up', 'Down']:
    sub = df[df['direction'] == grp]
    plt.scatter(sub['logFC'], sub['neglog10p'], s=28, c=colors[grp], label=grp, edgecolor='none')

# Axes labels and title
plt.xlabel('log2 fold change')
plt.ylabel('-log10 p-value')
plt.title('Volcanoplot - CES')

# Optional reference lines: p=0.05 and fc=0
plt.axhline(-np.log10(0.05), color='gray', lw=1, ls='--', alpha=0.6)
plt.axvline(0, color='gray', lw=1, ls='--', alpha=0.6)

# 6) Add labels only for top 5 features
for _, row in df.iterrows():
    name = str(row['Biomarkers'])
    if name in labels_top:
        # Offset text slightly to avoid overlap
        plt.text(row['logFC'] + 0.0002, row['neglog10p'] + 0.03, name, fontsize=8)

# Legend placement
plt.legend(title='Significance', frameon=False, loc='upper right')

plt.tight_layout()

out_png = 'volcano_CES.png'
plt.savefig(out_png, dpi=300)

out_png
