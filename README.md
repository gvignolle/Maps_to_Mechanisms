# Maps_to_Mechanisms

## Geospatial Factors Shape the Human Brain–Gut–Microbiome System — *Minimal Reproducible Analysis*

This repository ships a compact, execution‑ready subset of code to reproduce the manuscript’s core statistical findings: dose–response screens, cross‑modal integration (sPLS‑DA/DIABLO), and a publication‑style volcano plot. 

***

## What’s here

*   `minimal_code.R` — one‑file R workflow: data harmonization, univariate dose–response, nested CV for block sPLS‑DA (mixOmics), ROC/AUC, feature export, and SCFA↔enteroendocrine hub screen. Outputs include selected feature lists, hub CSV, and session info. 
*   `volcanoplot_python.py` — Python script to render a square volcano plot from a differential table (labels top‑5 features, saves `volcano_CES.png`). 

***

## Expected inputs (simple, rectangular CSV; one row per `SampleID`)

    data/
      labels.csv              # SampleID, CES_score, PBS_score, PopCS_score, CES_group, ZIP/Tract
      microbiome.csv          # SampleID, <species/taxa…>   (CLR/log CPM or VST)
      brain_fc.csv            # SampleID, <edge_1 … edge_P> (Fisher z)
      cytokines.csv           # SampleID, <GMCSF, IL17F, IL13, TNFA, VEGFA, …>
      brainstem_volumes.csv   # SampleID, <AAN_LC, AAN_PAG, AAN_VTA, AAN_PnO, …> (mm^3)
      scfa_hormones.csv       # SampleID, <acetate, propionate, butyrate, lactate, GLP1, PYY, insulin, Cpeptide, ghrelin>

These column names align with the minimal pipeline and manuscript terminology. 

***

## Quick start

### R (core analysis)

```r
# install once
install.packages(c("readr","dplyr","tidyr","pROC"))
# mixOmics from CRAN or Bioconductor
install.packages("mixOmics")
# run
Rscript minimal_code.R
```

Outputs (written to `outputs/`): (i) nested‑CV AUCs per block, (ii) final selected microbiome and brain‑edge panels, (iii) SCFA–hormone hub table, (iv) `sessionInfo.txt`. 

### Python (optional volcano plot)

```bash
pip install pandas numpy matplotlib openpyxl
python volcanoplot_python.py
```

Reads an Excel differential table and saves `volcano_CES.png` (top‑5 labels, FDR and FC guides). 

***

## Notes

*   This minimal code mirrors the manuscript’s analysis design (CalEnviroScreen strata, arousal‑network anatomy, cytokine gradients, and sPLS‑DA integration) using only the inputs above; full pipelines and data governance live with the primary study. 
*   Replace file paths as needed; figures may vary slightly by platform and random seeds. 

***

## Cite

**Vignolle GA *et al*.** “Geospatial Factors Shape the Human Brain–Gut–Microbiome System.” (see manuscript in repo for full author list and details). [Manuscript](http://something)

***

## Contact

For questions or data access, contact the corresponding authors listed in the manuscript. 

***

**License:** MIT
