# Differential Gene Expression Analysis of Hydrocortisone-Treated A549 Cells

### 🧬 Overview
This repository contains code, data, and results for identifying **differentially expressed genes (DEGs)** in **A549 cells treated with hydrocortisone** compared to **DMSO control** using **RNA-seq datasets from ENCODE**.

---

### 📂 Dataset
- **Control (DMSO):** ENCSR414IGI (5 replicates)
- **Hydrocortisone-treated:** ENCSR522CJR (5 replicates)

---

### ⚙️ Pipeline
1. Download quantification files
2. Merge counts by Ensembl Gene IDs
3. Perform normalization and DEG analysis using:
   - **DESeq2**
   - **edgeR**
4. Generate visualizations:
   - Volcano plots
   - MA plots
   - Heatmaps
5. Identify overlapping DEGs between methods

---

### 🧠 Results Summary
| Metric | DESeq2 | edgeR | Overlap |
|--------|--------|--------|----------|
| Significant DEGs | 715 | 632 | 614 |
| Up-regulated | 553 | 487 | 475 |
| Down-regulated | 162 | 145 | 139 |
| Overlap (%) | **91.17%** |  |  |

---

### 🔍 Key Findings
- Strong induction of **FKBP5**, **TSC22D3 (GILZ)**, **ANGPTL4**, and **PER1**
- Repression of **CCL2**, **IL11**, **NRG1**, and **WNT7B**
- Results are consistent between DESeq2 and edgeR
  
---

### 🧾 Citation
If you use this repository, please cite:
> Arif, A. (2025). *Differential Gene Expression Analysis of Hydrocortisone-Treate*
