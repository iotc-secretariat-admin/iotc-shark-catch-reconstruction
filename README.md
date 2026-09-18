# IOTC Blue Shark and Shortfin Mako Historical Catch Reconstruction

This repository contains the analytical workflow used to reconstruct historical retained catches of blue shark (*Prionace glauca*, BSH) and shortfin mako (*Isurus oxyrinchus*, SMA) in the Indian Ocean for IOTC-WPEB22.

The analysis accompanies:

> Rice, J. 2026.  Reconstructing Historical Catch Series for Blue Shark (*Prionace glauca*) and Shortfin Mako (*Isurus oxyrinchus*) in the Indian Ocean: Methodological Approaches and Uncertainty Quantification.  IOTC-2026-WPEB22(AS)-30.

## Overview

Historical shark catch records contain several distinct forms of missing or taxonomically unresolved information. The reconstruction therefore does not treat all missing species-specific catch as a single estimation problem.

The analysis distinguishes among:

1. **Species-specific reported catch**
2. **Taxonomically unresolved shark catch**
3. **Completely unreported shark catch**

These components require different estimation procedures. Taxonomically unresolved catches are reconstructed using hierarchical species-composition information, whereas complete non-reporting is addressed using donor-based species-to-target catch relationships. Alternative reconstruction scenarios and statistical sensitivity analyses are used to characterize uncertainty associated with these historical allocations.

A central objective of the repository is to make the complete analytical workflow **reproducible and auditable**. Consistent with the project Terms of Reference, the repository is structured so that data collation, processing, reconstruction, sensitivity analyses, and scenario construction can be traced from their original inputs through intermediate analytical steps to the final reconstructed catch series. Assumptions and methodological choices are implemented explicitly in the analysis code rather than being applied manually to final results.

The repository therefore serves both as the computational basis for the WPEB22 analysis and as a documented record of how the reconstructed catch histories were produced. Because IOTC source datasets may be revised over time, retaining the source-data version, analytical code, and scenario definitions together is important for distinguishing changes in methodology from changes in the underlying data.

## Using this repository

The repository is organized as a reproducible analysis rather than simply as a collection of scripts. Users interested only in the results can review the working paper and final outputs. Users wishing to examine or reproduce the analysis should retain the repository directory structure and run the workflow from the repository root.

For users unfamiliar with Git or GitHub, a **repository** is a version-controlled project containing the analysis code, documentation, and the history of changes to those files. GitHub provides the online copy of the repository, while Git can be used to create and maintain a local copy on your computer.

The recommended workflow is:

1. **Obtain the repository.** Users familiar with Git can clone the repository. Alternatively, GitHub provides an option to download a snapshot of the repository as a ZIP file. Cloning is preferred when the analysis will be updated or changes need to be tracked.

2. **Retain the directory structure.** Scripts use repository-relative paths such as `data/raw/`, `R/`, and `outputs/`. Moving individual scripts or data files outside this structure may prevent the workflow from running correctly.

3. **Place or verify the source data.** The analysis uses the IOTC retained-catch (RC), geo-referenced catch (CA), and effort (EF) datasets identified in the repository documentation. Source-data filenames and versions should be checked before reproducing the analysis because subsequent IOTC data revisions may change the resulting reconstruction.

> > **IMPORTANT:** The IOTC retained-catch (RC) dataset used in this analysis is included in `data/raw/rc/`. The catch-and-effort datasets ([CA and EF](https://iotc.org/sites/default/files/documents/2026/08/IOTC-DATASETS-2026-08-13-CE-1952-2025.zip)) are not version controlled because of their size and must be obtained separately from IOTC. The analysis used the 13 August 2026 releases covering 1952–2025. After downloading and unzipping the archive, the CA and EF files must be placed at the following locations:
>
> - `data/raw/ca/IOTC-DATASETS-2026-08-13-CA-1952-2025.csv`
> - `data/raw/ef/IOTC-DATASETS-2026-08-13-EF-1952-2025.csv`
>
> The production workflow expects these filenames and directory locations.

4. **Run the production workflow.** From R with the working directory set to the repository root, run:

   ```r
   source("R/run_reconstruction.R")
   ```

   This script executes the principal data-processing and reconstruction components in their required order and produces the final reconstructed catch outputs.  
   
   **Note**: Running the complete reconstruction may take considerable time, depending on the workstation. The statistical sensitivity analyses, particularly repeated GAM cross-validation, account for much of the computational time.  

5. **Review the outputs and diagnostics.** Final reconstructed catch series are written to `outputs/`, while species-specific and diagnostic results are retained in the corresponding subdirectories. Supporting exploratory and validation analyses are contained separately under `R/exploration/` and `R/diagnostics/`.

The `main` branch should be treated as the stable version of the analysis. Users modifying or extending the workflow are encouraged to make changes on a separate branch, use small and clearly described commits, and merge reviewed changes back into `main`. This preserves a transparent history of analytical changes and helps distinguish modifications to data processing, assumptions, reconstruction methods, and presentation outputs.

For readers new to GitHub, the official GitHub documentation provides an introductory walkthrough of repositories, branches, commits, and pull requests. More detailed recommendations for maintaining repositories are available in GitHub's repository best-practices guidance.

**GitHub resources:**

* GitHub introductory tutorial: [Hello World](https://docs.github.com/en/get-started/using-github/hello-world)
* GitHub repository guidance: [Best practices for repositories](https://docs.github.com/en/repositories/creating-and-managing-repositories/best-practices-for-repositories)

## Reproducibility and traceability

The intended analytical chain is:

```text
IOTC source data
      ↓
data import and code harmonization
      ↓
reconstruction strata
      ↓
reporting-state and species-composition analyses
      ↓
BSH and SMA reconstruction components
      ↓
alternative scenarios and sensitivity analyses
      ↓
diagnostics and validation
      ↓
final annual reconstructed catch series
```

Where practical, derived quantities are generated by the scripts rather than edited manually. This allows users to trace a reported result back through the reconstruction procedure to the relevant input data, donor hierarchy, assumptions, and scenario definition.

The primary reproducibility test is to start a clean R session from the repository root and run:

```r
rm(list = ls())
gc()

source("R/run_reconstruction.R")
```

A successful clean-session run helps ensure that the reconstruction does not depend on objects created interactively or retained from an earlier R session.

## Author and contact

**Joel Rice**
Rice Marine Analytics
Fisheries stock assessment and quantitative fisheries science
Email: **ricemarineanalytics@gmail.com**

Questions regarding the reconstruction methodology, analytical assumptions, or implementation of the workflow are welcome. When reporting an issue, please identify the relevant script, data version, and—where applicable—the error message or output being examined.

## Repository structure 

The repository is organized to separate source data, data processing, reconstruction, diagnostics, exploratory analyses, outputs, and documentation. This structure is intended to make the analytical workflow easier to navigate and to preserve a clear connection between source inputs, analytical code, intermediate products, and final reconstructed catch series.

Users reproducing the analysis should retain this directory structure because the production scripts use repository-relative paths. The R/reconstruction/ directory contains the core BSH and SMA reconstruction procedures, while R/diagnostics/ and R/exploration/ contain supporting analyses that are useful for evaluating the reconstruction but are not necessarily required to generate the final catch series. The principal assessment-ready results are written to outputs/, and the working-paper source and supporting documentation are maintained in docs/.

 

```text
.
├── README.md
│
├── data/
│   └── raw/
│       ├── ca/
│       │   └── IOTC-DATASETS-2026-08-13-CA-1952-2025.csv
│       ├── ef/
│       │   └── IOTC-DATASETS-2026-08-13-EF-1952-2025.csv
│       ├── rc/
│          └── IOTC-DATASETS-2026-08-13-RC-RAW-1950-2025.csv
│        
│
├── R/
│   ├── processing/
│   │   ├── Initial_data_import.R
│   │   ├── harmonize_iotc_codes.R
│   │   └── build_reconstruction_strata.R
│   │
│   ├── diagnostics/
│   │   ├── compare_rc_ca_shark_composition.R
│   │   ├── compare_rc_ca_sma_composition.R
│   │   └── plot_bsh_rc_history.R
│   │
│   ├── exploration/
│   │   ├── data inventory and reporting diagnostics
│   │   ├── spatial reporting maps
│   │   ├── spatial catch and effort analyses
│   │   └── SMA data inventory
│   │
│   ├── reconstruction/
│   │   ├── bsh_shark_nei_disaggregation.R
│   │   ├── bsh_nei_scenarios.R
│   │   ├── bsh_nei_historical_sensitivity.R
│   │   ├── bsh_ratio_reconstruction.R
│   │   ├── bsh_gam_species_composition.R
│   │   ├── reconstruct_bsh_scenarios.R
│   │   ├── sma_reporting_states.R
│   │   ├── sma_hscr_species_composition.R
│   │   ├── sma_hscr_nonreporting.r
│   │   ├── sma_gam_species_composition.R
│   │   └── reconstruct_sma_scenarios.R
│   │
│   └── utils/
│       └── check_repo.R
│
├── outputs/
│   ├── bsh/
│   │   └── gam/
│   │       ├── summary and comparison tables
│   │       └── figures/
│   │
│   ├── sma/
│   │   ├── reconstruction summaries and intermediate objects
│   │   └── gam/
│   │       ├── model diagnostics and comparison tables
│   │       └── figures/
│   │
│   ├── spatial/
│   │   ├── coverage summaries
│   │   ├── fishing-ground lookup files
│   │   └── figures/
│   │
│   └── final_BSH_SMA_annual_reconstructed_catch.csv
│
└── docs/
    ├── Reconstructing_hist_BSH_SMA_catches_IOTC_updated.qmd
    ├── refs.bib
    ├── ices-journal-of-marine-science.csl
    ├── journal_doc_template_ENG.docx
    └── rendered DOCX  manuscript file 
     
```

## Main analytical workflow

The broad sequence is:
 

```text
Raw IOTC RC / CA / EF data
        |
        v
R/processing/
  import, harmonize codes, construct reconstruction strata
        |
        +---------------------------+
        |                           |
        v                           v
R/diagnostics/                R/exploration/
  RC–CA comparisons             reporting history
  composition checks            spatial patterns
  historical checks             effort overlap
        |                           |
        +-------------+-------------+
                      |
                      v
             R/reconstruction/
                      |
          +-----------+-----------+
          |                       |
          v                       v
        BSH                       SMA
          |                       |
          |              reporting-state classification
          |              MAK allocation
          |              SKH HSCR allocation
          |              non-reporting reconstruction
          |              GAM sensitivities
          |                       |
          +-----------+-----------+
                      |
                      v
                 outputs/
                      |
                      v
 final_BSH_SMA_annual_reconstructed_catch.csv
                      |
                      v
                    docs/
               Quarto working paper
```

## Reconstruction production pipeline

```text
1. R/processing/Initial_data_import.R
        |
        v
2. R/processing/harmonize_iotc_codes.R
        |
        v
3. R/processing/build_reconstruction_strata.R
        |
        +-----------------------------------+
        |                                   |
        v                                   v
      BSH branch                          SMA branch
        |                                   |
4a. bsh_shark_nei_disaggregation.R       4b. sma_reporting_states.R
        |                                   |
5a. bsh_nei_scenarios.R                  5b. sma_hscr_species_composition.R
        |                                   |
6a. bsh_nei_historical_sensitivity.R     6b. sma_hscr_nonreporting.R
        |                                   |
7a. bsh_ratio_reconstruction.R           7b. sma_gam_species_composition.R
        |                                   |
8a. bsh_gam_species_composition.R        8b. reconstruct_sma_scenarios.R
        |
9a. reconstruct_bsh_scenarios.R
        |
        +------------------+----------------+
                           |
                           v
             final_BSH_SMA_annual_reconstructed_catch.csv

```

### Diagnostics are not necessarily required to run the reconstruction

```text
R/diagnostics/
    compare_rc_ca_shark_composition.R
    compare_rc_ca_sma_composition.R
    plot_bsh_rc_history.R

R/exploration/
```


## Reconstruction framework

### Taxonomic aggregation

Catch reported under unresolved shark categories already contains
information on total shark catch. The remaining problem is allocation
among species.

Generic Sharks NEI (`SKH`) are allocated using **Hierarchical
Species-Composition Reconstruction (HSCR)**.

For donor set $d(i)$, the focal-species proportion is

```math
\hat{p}_{s,d(i)}
=
\frac{
\sum_{j \in d(i)} C_{s,j}
}{
\sum_{j \in d(i)} C_{\mathrm{identified\ shark},j}
}.
```

Generic shark catch is then allocated as

```math
\hat{C}_{s,i}^{\mathrm{SKH}}
=
C_{\mathrm{SKH},i}
\hat{p}_{s,d(i)}.
```

Two principal HSCR scenarios are retained:

- **N1** — more restrictive empirical donor hierarchy
- **N2** — broader empirical donor support


### Shortfin mako: unresolved mako catch

SMA requires an additional taxonomic reconstruction because catches may
be reported as mako sharks NEI (`MAK`).

For donor set \(d(i)\), the SMA proportion of identified mako catch is

```math
\hat{p}_{\mathrm{SMA},d(i)}
=
\frac{
\sum_{j \in d(i)} C_{\mathrm{SMA},j}
}{
\sum_{j \in d(i)}
\left(
C_{\mathrm{SMA},j}
+
C_{\mathrm{LMA},j}
\right)
}.
```

Unresolved mako catch is then allocated as

```math
\hat{C}_{\mathrm{SMA},i}^{\mathrm{MAK}}
=
C_{\mathrm{MAK},i}
\hat{p}_{\mathrm{SMA},d(i)}.
```


### Complete non-reporting

Where target catch is reported but shark catch is absent entirely,
species-specific catch is reconstructed using empirically supported
shark-to-target catch relationships.

For donor set \(d(i)\), the species-to-target catch ratio is

```math
\hat{R}_{s,d(i)}
=
\frac{
\sum_{j \in d(i)} C_{s,j}
}{
\sum_{j \in d(i)} C_{\mathrm{target},j}
}.
```

Catch in recipient stratum \(i\) is then reconstructed as

```math
\hat{C}_{s,i}
=
C_{\mathrm{target},i}
\hat{R}_{s,d(i)}.
```

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

```math
\mathrm{logit}(E[p_s])
=
\mathrm{Fishery}
+
s(\mathrm{Year}, by=\mathrm{Fishery})
+
s(\mathrm{Fleet}, RE)
+
\sum_k s(\mathrm{Target\ composition}_k).
```

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



 



## Blue shark reconstruction

The final BSH reconstruction contains:

- reported BSH;
- BSH allocated from SKH; and
- BSH estimated for complete non-reporting.

HSCR N1 and N2 provide the principal empirical reconstruction range.

A beta-GAM species-composition model is retained as a sensitivity.

## Shortfin mako reconstruction

The final SMA reconstruction contains:

- reported SMA;
- empirical allocation of MAK between SMA and LMA;
- HSCR allocation of SKH; and
- supported reconstruction of complete non-reporting.

The SMA GAM analyses are retained as sensitivities or diagnostics rather
than preferred reconstructions.

## Principal output

The main assessment-ready output is: 

> outputs/final_BSH_SMA_annual_reconstructed_catch.csv

It contains annual:

- reported BSH;
- BSH HSCR N1;
- BSH HSCR N2;
- recommended BSH;
- reported SMA;
- SMA HSCR N1;
- SMA HSCR N2; and
- recommended SMA.

The recommended annual estimate is currently reported as the midpoint
between N1 and N2, while both empirical reconstruction scenarios are
retained explicitly.

## Reproducibility

Historical IOTC catch records are revised through time. Results therefore
depend on both:

+ the reconstruction method; and
+ the version of the underlying IOTC data.

Assessment-specific input datasets and analytical outputs should be
retained so changes caused by source-data revisions can be distinguished
from changes in reconstruction methodology.


## Manuscript

The working paper source is:

> docs/Reconstructing_hist_BSH_SMA_catches_IOTC_updated.qmd

The bibliography is:

> docs/refs.bib

Rendered  DOCX file is excluded from version control and can be
regenerated from the Quarto source. Will need minor formatting prior to submission. For example the 'cover sheet' questions are produced and answered in the Appendix.
