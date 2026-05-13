# SPRUCE Litter Decomposition Project 2015 - 2018 

This repository contains the code for the analyses and figures presented in the manuscript:
 
**Litter Identity Drives Molecular Trajectories of Peatland Organic Matter Decomposition Under Warming and Elevated CO<sub>2</sub>**
 
Data are available in this GitHub repository and in the OSF repository (DOI: *to be added upon publication*). 

---
 
## Citation
 
*To be added upon publication.*

---
## Overview
 
Litter decomposition controls the transfer of plant-derived carbon into soil organic matter, yet the molecular mechanisms driving decomposition remain poorly understood, particularly across diverse litter types under changing climate conditions. This study examines how litter identity, warming, and elevated CO<sub>2</sub> shape molecular-level organic matter transformations in a peatland ecosystem.
 
Six litter types representing dominant plant functional groups of the S1 bog were incubated at the [SPRUCE](https://mnspruce.ornl.gov/) (Spruce and Peatland Responses Under Changing Environments) whole-ecosystem experiment in northern Minnesota across three warming levels (+0, +4.5, +9°C) and two CO<sub>2</sub> concentrations (ambient and +500 ppm), with collections at 0, 0.5, 1, and 2 years. Litter chemistry was characterized using bulk elemental analysis, FTIR spectroscopy, and Fourier Transform Ion Cyclotron Resonance Mass Spectrometry (FT-ICR MS), processed via the [MetaboDirect](https://github.com/Coayala/MetaboDirect) pipeline. 
 
**Litter types (ordered by increasing decomposition rate):** MAG (*Sphagnum magellanicum*), SPR (spruce roots), LTR (Labrador tea roots), ANG (*Sphagnum angustifolium*), SPL (spruce needles), LTL (Labrador tea leaves).

-----
## Scripts
 
### Bulk Chemistry
- **`Fig1_SupFig1-2-BulkCharacterization.R`** — Decomposition rates (k), initial C:N, %C, %N, %P, and FTIR-derived carbohydrate and aromatic content across litter types and timepoints.
- **`SupFig3-BulkvsFTICR.R`** — Comparison of bulk FTIR vs. FT-ICR MS carbohydrate and aromatic content.
### Molecular Composition & Trajectories
- **`Fig2-Overall_Change_Trajectory.R`** — NMDS ordination and PCoA-based trajectory analysis (trajectory length and speed) using `ecotraj`.
- **`Fig3_transformations.R`** — Pairwise molecular transformation profiles over time following Fudyma et al. (2019).
- **`Fig4_Permanova_Table1.R`** — PERMANOVA of molecular class composition by litter, time, temperature, and CO<sub>2</sub>; per-litter bubble plot.
- **`SupFig4-5-6_OverallClasses_Composition_NOSC_GFE.R`** — Overall compound class shifts, elemental composition, and thermodynamic indices (NOSC, GFE) over time across all litter types.
### Molecular Diversity & Turnover
- **`Fig5_ANG_ShannonDiversity.R`** — Shannon diversity of compound classes over time for *S. angustifolium*.
- **`Fig6_VK_ANG_SupFig9.R`** — Van Krevelen turnover diagrams and thermodynamic index trajectories for *S. angustifolium*.
- **`Fig7_Shannon_Litter_Comparison.R`** — Shannon diversity trajectories across all litter types and timepoints.
- **`SupFig7_MAG_ShannonDiversity.R`** — Shannon diversity of compound classes over time for *S. magellanicum*.
- **`SupFig8-10_VK_MAG_NOSC_GFE.R`** — Van Krevelen turnover and thermodynamic indices for *S. magellanicum*.
- **`SupFigs11to22_Temp&CO2 Effect on All Litters.R`** — Shannon diversity broken down by temperature and CO<sub>2</sub> treatment for all six litter types.
### Predictors of Decomposition
- **`Fig8_CorrelationWithMass.R`** — Spearman correlations between compound-class changes and mass loss, across and within litter types.
- **`SupTable2_Index.R`** — Linear regressions of T0 molecular indices (AI, AI_mod, GFE, NOSC, DBE, DBE-O, C:N) against decomposition rate (k), with ranked R² output.

 
