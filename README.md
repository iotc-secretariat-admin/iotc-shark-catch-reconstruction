# IOTC Blue Shark and Shortfin Mako Historical Catch Reconstruction

This repository contains the analytical workflow used to reconstruct
historical retained catches of blue shark (*Prionace glauca*, BSH) and
shortfin mako (*Isurus oxyrinchus*, SMA) in the Indian Ocean for
IOTC-WPEB22.

The analysis accompanies:

> Rice, J. 2026.  Reconstructing Historical Catch Series for Blue Shark
> (Prionace glauca) and Shortfin Mako (Isurus oxyrinchus) in the Indian
> Ocean: Methodological Approaches and Uncertainty Quantification. 
> IOTC-2026-WPEB22(AS)-30.


## Overview

Historical shark catch records contain several different forms of missing
or unresolved information. The reconstruction therefore does not treat
all missing species-specific catch as a single problem.

The analysis distinguishes between:

1. **Species-specific reported catch**
2. **Taxonomically unresolved shark catch**
3. **Completely unreported shark catch**

These require different estimation procedures.


## Reconstruction framework

### Taxonomic aggregation

Catch reported under unresolved shark categories already contains
information on total shark catch. The remaining problem is allocation
among species.

Generic Sharks NEI (`SKH`) are allocated using **Hierarchical
Species-Composition Reconstruction (HSCR)**.

For donor set \(d(i)\),

$$
\hat{p}_{s,d(i)}
=
\frac{
  \sum_{j \in d(i)} C_{s,j}
}{
  \sum_{j \in d(i)} C_{\mathrm{identified\ shark},j}
}.
$$

Generic shark catch is then allocated as

\[
\hat{C}_{s,i}^{SKH}
=
C_{SKH,i}
\hat{p}_{s,d(i)}.
\]

Two principal HSCR scenarios are retained:

- **N1** — more restrictive empirical donor hierarchy
- **N2** — broader empirical donor support


### Shortfin mako: unresolved mako catch

SMA requires an additional taxonomic reconstruction because catches may
be reported as mako sharks NEI (`MAK`).

The SMA proportion of identified mako catch is estimated as

\[
\hat{p}_{SMA,d(i)}
=
\frac{
  \sum C_{SMA}
}{
  \sum (C_{SMA} + C_{LMA})
}.
\]

This proportion is used to allocate unresolved `MAK` catch between
shortfin and longfin mako.


### Complete non-reporting

Where target catch is reported but shark catch is absent entirely,
species-specific catch is reconstructed using empirically supported
shark-to-target catch relationships.

For donor group \(d(i)\),

$$
\hat{R}_{s,d(i)} 
=
\frac{ \sum C_s}{
  \sum C_{\mathrm{target}}
}.
$$

Catch in recipient stratum $(i)$ is then

\[
\hat{C}_{s,i}
=
C_{\mathrm{target},i}
\hat{R}_{s,d(i)}.
\]

Donor selection is hierarchical and retains information on temporal
support:

- within donor period;
- backcast; and
- forward extrapolation.

Strata lacking a suitable donor are retained as unresolved rather than
being reconstructed using increasingly dissimilar fisheries.


## Statistical sensitivities

Generalized additive models (GAMs) were evaluated as statistical
alternatives to empirical HSCR allocation.

The response is species composition within identified shark catch and is
modelled using beta regression.

The general model structure is

\[
\mathrm{logit}(E[p_s])
=
\mathrm{Fishery}
+
s(\mathrm{Year}, by=\mathrm{Fishery})
+
s(\mathrm{Fleet}, RE)
+
\sum_k s(\mathrm{Target\ composition}_k).
\]

Model support is evaluated primarily through:

- leave-one-fleet-out cross-validation;
- MAE;
- RMSE;
- bias;
- Spearman rank correlation;
- fishery-specific predictive performance;
- temporal extrapolation; and
- comparison with empirical HSCR estimates.

The GAM analyses are treated as **sensitivities or diagnostics**, not
automatically as preferred reconstructions.


## Uncertainty

Uncertainty is characterized using several complementary approaches:

- HSCR N1 versus N2 donor hierarchies;
- GAM versus empirical HSCR reconstruction;
- leave-one-fleet-out prediction;
- temporal-support classification;
- fishery-specific diagnostics; and
- explicit retention of unsupported strata as unresolved.

These represent structural and empirical uncertainty rather than a
single probabilistic confidence interval.


## Data sources

Three principal IOTC datasets are used:

- **RC** — retained catch
- **CA** — geo-referenced catch and effort
- **EF** — effort

The primary reconstruction unit is:

`Fleet × Fishery × Year`

Spatial CA and EF data are used primarily as diagnostics of reporting
history, spatial coverage, and fishing activity rather than to directly
allocate reconstructed catch.


## Repository structure

```text
.
├── data/
│   └── Input data or local data links
│
├── metadata/
│   └── Data dictionaries and supporting metadata
│
├── R/
│   ├── functions/
│   │   └── Shared functions
│   │
│   ├── import/
│   │   └── Data import and standardization
│   │
│   ├── processing/
│   │   └── Construction of common analytical strata
│   │
│   ├── reconstruction/
│   │   ├── BSH reconstruction
│   │   ├── SMA reconstruction
│   │   ├── HSCR allocation
│   │   ├── complete non-reporting reconstruction
│   │   └── GAM sensitivity analyses
│   │
│   ├── diagnostics/
│   │   └── Cross-validation and reconstruction diagnostics
│   │
│   └── plotting/
│       └── Manuscript and diagnostic figures
│
├── outputs/
│   ├── bsh/
│   ├── sma/
│   ├── spatial/
│   └── final_BSH_SMA_annual_reconstructed_catch.csv
│
└── docs/
    ├── Reconstructing_hist_BSH_SMA_catches_IOTC_updated.qmd
    └── refs.bib

```