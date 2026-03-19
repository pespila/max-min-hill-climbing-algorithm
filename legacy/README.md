# Legacy R/C++ Implementation

This directory contains the **original R package implementation** of the Max-Min Hill-Climbing (MMHC) algorithm, preserved for reference.

The code was written by Michael Bauer (2014) as an R package using Rcpp for C++ acceleration. It has since been superseded by the Python implementation in the `mmhc/` directory.

## Structure

- **`R/`** — R source code (`data.R`, `mmhc.R`)
- **`src/`** — C++ source code via Rcpp (`mmhc.cpp`, `mmhc.h`, `mmpc.cpp`, `rcpp_module.cpp`)
- **`man/`** — Documentation and the original report (`mmhc.Rd`, `report.pdf`)
- **`tests/`** — Test script (`test.R`)
- **`DESCRIPTION`** / **`NAMESPACE`** — R package metadata

## Source

Extracted from commit `8a43f8b`, the last commit before the Python migration.
