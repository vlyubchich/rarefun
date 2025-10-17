# rarefun: Functions for Rare Events Analysis

[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/rarefun)](https://CRAN.R-project.org/package=rarefun)
[![R-CMD-check](https://github.com/vlyubchich/rarefun/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/vlyubchich/rarefun/actions/workflows/check-standard.yaml)

## Overview

`rarefun` is an R package that provides a collection of functions for detecting and analyzing rare events in datasets. The package includes implementations of various anomaly detection methods and tools for spatial and temporal analysis of rare events.

## Features

- **Anomaly Detection Methods:**
  - Isolation Forest (`rare_iforest`)
  - DBSCAN clustering (`rare_dbscan`)
  - Residual-based analysis (`rare_residuals`)

- **Spatial Analysis:**
  - Spatial point matching (`match_spatial_points`)
  - Geographic data processing

- **Time Series Tools:**
  - Time series extraction (`extract_ts`)
  - Temporal pattern analysis

- **Statistical Methods:**
  - Contingency table chi-square tests (`contingency_chisq`)
  - Goodness-of-fit tests for quantile regression (`gof_qr`)
  - Matthews Correlation Coefficient (`mcc`)
  - Partial quantile random forest (`partial_qrf`)
  - Kneedle algorithm for knee point detection (`kneedle`)

## Installation

You can install the development version of rarefun from GitHub:

```r
# install.packages("devtools")
devtools::install_github("vlyubchich/rarefun")
```

## Usage

```r
library(rarefun)

# Example: Detect anomalies using Isolation Forest
set.seed(123)
data <- matrix(rnorm(1000), nrow = 500)
data[1:5, ] <- data[1:5, ] + 10  # Add some outliers

result <- rare_iforest(data, ntrees = 100, threshold = 0.7)
table(result$is_anomaly)
```

## Authors

- **Vyacheslav Lyubchich** (Maintainer) - [ORCID: 0000-0001-7936-4285](https://orcid.org/0000-0001-7936-4285)
- **Geneviève Nesslage** - [ORCID: 0000-0003-1770-6803](https://orcid.org/0000-0003-1770-6803)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Bug reports

If you encounter any bugs or issues, please report them at: https://github.com/vlyubchich/rarefun/issues
