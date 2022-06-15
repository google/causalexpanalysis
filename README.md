# CausalExpAnalysis R Package

#### Causal and Experiment Analysis

The CausalExpAnalysis R package contains functions to perform causal and
experiment analysis, i.e., analysis for the causal inference of both experiments
and observational studies, within a single framework, with a unifying theme of
evaluating covariate balance. Three kinds of samples can be created—unadjusted,
weighted, and matched—and two kinds of estimators of treatment effects—an
analysis of variance (ANOVA) estimator, i.e., a simple difference-in-means
estimator, and an analysis of covariance (ANCOVA) estimator, i.e., a regression
adjustment estimator, with either equal or unequal slopes assumption for the
covariates. Contains helper functions for statistical inference and modeling,
visualization of results, and data manipulation.

#### Installation

```r
install.packages("devtools")
library(devtools)
devtools::install_github("google/CausalExpAnalysis")
library(CausalExpAnalysis)
```
