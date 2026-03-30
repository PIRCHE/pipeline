# 2026 ASHI Clinical Informatics Workshop: Connect to the PIRCHE API with R

Materials for automating Solid Organ Transplant (SOT) PIRCHE-T2 and PIRCHE-B score calculations using R and the PIRCHE REST API. 



### Folder Contents

* **`ACIW_PIRCHE-intro.pdf`** Introductory slide deck serving as a beginner's primer. Covers the basics of SOT PIRCHE scores, API POST requests, and JSON formatting before transitioning to the code.

* **`test_data.csv`** Sample input dataset of mock clinical cases.

* **`pirche-sot-wkshp.Rmd`** Interactive R Markdown notebook. Provides a step-by-step tutorial of the API workflow, detailed breakdowns of key functions, and side-by-side comparisons with standard Bash `cURL` requests.

* **`pirche-sot-score-requests.R`** Consolidated R script. Executes data preparation, API payload construction, and score retrieval. Includes a function to construct a deep link that opens the TxPredictor web application with pre-filled case data.

---

### Prerequisites

Running these scripts requires R and the following packages:
```
install.packages(c("tidyverse", "httr2", "immunogenetr", "curl"))
```

### Usage

1. Clone the parent repository or download this folder locally.
2. Obtain your personal PIRCHE API credentials. Contact `support@pirche.com` for assistance.
3. Open `pirche-sot-wkshp.Rmd` (for the tutorial) or `pirche-sot-score-requests.R` (for the automated pipeline) in RStudio.
4. Input your credentials where prompted and execute the code chunks sequentially.

---

### License and Citation

This code is open-source under the MIT License (see the root repository `LICENSE` file).

If you use or adapt these materials in your research or clinical pipeline, please cite:

> <small>**Mehler, Hilary.** (2026). *Connect to the PIRCHE API with R.* Materials presented at the 2026 ASHI Clinical Informatics Workshop. GitHub repository: https://github.com/PIRCHE/pipeline</small>

