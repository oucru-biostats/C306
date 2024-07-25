> [!IMPORTANT]
> THIS BRANCH IS FEATURE-FROZEN, AIMING FOR REPRODUCIBILITY IN 27TB ANALYSIS. ONLY BUG FIXES SHOULD BE ACCEPTED


# C306

This package provides a collection of useful R codes for data analysis
from the Biostatistics group at OUCRU-HCM.

## Installation

You can install C306 from github with:

``` r
# install.packages("devtools")
devtools::install_github("oucru-biostats/C306@feature-rbind")
```

## List of functions

### Main functions

These functions’ purpose are to generate report tables, outputting a
flextable or a intermediate object of class `ss_table` for post-hoc
manipulation.

- `sstable.baseline` function to generate baseline table 1.
- `sstable.survcomp` function to generate dual arm descriptive and
  inferential analysis for time-to-event outcome, supporting *Cox pH*
  and *Restricted mean surival time* models
- `sstable.survcomp.subgroup` extending `sstable.survcomp` for subgroup
  analyses
- `sstable.ae` function to generate adverse event tables with $\chi^2$
  or Fisher test for comparision.

### sstable manipulation

These function provides helper to create ss_table-esque output, with
supporting function to convert them to flextable or huxtable.

- `ssformat` function to convert an arbitrary matrix to an sstable
- `ss_flextable` function to convert an sstable to a flextable
- `ss_huxtable` function to convert an sstable to a huxtable
- `ft_sstheme` function to decorate an arbitrary flextable to have an
  sstable-esque theme
- `ht_theme_markdown` and `ht_theme_kable` functions to decorate
  huxtable to follow a sstable-esque theme
- `as_sstable` converting objects to sstable

## Plotting functions

- `gg_boxcox` same with `MASS::boxcox` but in ggplot
- `ggsurvfit2` same as `ggsurvfit::ggsurvfit` but with strata separated.
  Note that it relies on `tidy_survfit2` and does not work with
  `ggsurvfit::survfit2` at the moment.
- `gg_ajsurvplot` plots one Aalen-Johansen curve for competing risks for
  the event of interest. Note that this is based on `surviminer` due to
  its flexibility
- `gg_ajsurvplot2` plots **two** Aalen-Johansen curves for main risk
  (from bottom) and competing risk (from top). Note that this is based
  `ggsurvfit` with some limitations. However, it has the ability to
  return a dataset for more flexibility in plotting with ggplot.

### OUCRU function written by Lam PK

- `import.info` and `convert.info` import and convert OUCRU dictionary
  to C306 style dict to use for data inspection
- `insepct.data` function to check data for error based on a dictionary
  by `convert.info`
- `myformat.data` function to reformat the data to readable form base on
  the dictionary by `convert.info`
- `labAE` function to determine a laboratory AE based on a threshold

### Helper function written by Trinh Dong

- `logist_summary` based on a function by OUCRU Biostats group,
  reporting the OR of logistic regression
- `subgroup_effect` extends `logist_summary` to report the subgroup OR
  of variables that have interaction with a common covariate
  (i.e. treatment arm).
- `mutate_f` and `summarise_f` perform multiple assignment for
  data.frame by `c(a,b):=list(v_a, v_b)`
- `simple_relevel` relevels a factor correposnding to the level of
  another factor
- `simple_recode` pattern-based recoding of data by in-place replacement
- `tidy_survfit2` same as but with strata saved for easier plotting.
  Note that it, at the moment, does not work with `ggsurvfit::survfit2`
  ironically.
