"""# Project Documentation: QGenomix Phase 1

**Work Period Covered:** 26 May 2026  
**Project Root:** `D:\\QGenomix_mem1_step1`  
**Objective:** Automate the download, quality control, and linkage disequilibrium (LD) pruning of 1,000 Genomes autosomal data (Chromosomes 1–22), and extract a top-1,000 Type 2 Diabetes (T2D) risk-associated SNP matrix for downstream analysis.

---

## 1. Environment & Tooling
* **Tool:** PLINK v2.0.0-a.7.1 64-bit (4 May 2026)
* **Location:** `\\tools\\plink2`
* **Purpose:** Genomic data reading, filtering, pruning, merging, and exportation.

---

## 2. Data Pipeline & Quality Control (QC)

### Automation Script
* **Path:** `\\run_download_qc_prune.ps1`
* **Scope:** Chromosomes 1 to 22
* **Output Directory:** `\\data\\interim\\qc_by_chr`

### Per-Chromosome Processing Steps
1. Download compressed genotype (`.pgen.zst`) and no-annotation variant (`.pvar.zst`) files.
2. Decompress `.pgen.zst` and link the corrected sample file (`.psam`).
3. Apply rigorous QC filters and resolve duplicate variant IDs.
4. Perform LD pruning.
5. Clear temporary decompressed files to optimize disk space.

### QC and Pruning Thresholds
| Parameter | Filter Value |
| :--- | :--- |
| **Minor Allele Frequency (MAF)** | &ge; 0.01 |
| **Hardy-Weinberg Equilibrium (HWE)** | $p\\text{-value} \\ge 1\\times 10^{-6}$ |
| **Variant Missingness** | &le; 0.05 |
| **Sample Missingness** | &le; 0.1 |
| **Allele Type Constraints** | Biallelic, A/C/G/T SNPs only |
| **LD Pruning Window** | 50 variant window, 5 variant step, $r^2 > 0.2$ |

---

## 3. Key Issues Resolved

### Duplicate Variant IDs
During the initial processing run, PLINK flagged errors regarding non-unique variant IDs. 
* **Resolution:** Updated the automation script to dynamically construct unique synthetic IDs using the pattern `chromosome:position:reference_allele:alternate_allele` (e.g., `1:123456:A:G`) and dropped unresolved duplicates.
* **PLINK2 Flags Added:**
  ```bash
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 missing --rm-dup exclude-all
  
---

## 4. Final Datasets & Output Summary

### Merged Master Dataset
Individual cleaned chromosome datasets were merged into a single autosomal master dataset[cite: 1].
* **Location:** `\data\interim\merged\all_autosomes_pruned.*` (`.pgen`, `.pvar.zst`, `.psam`)[cite: 1]

> **Dataset Demographics:**
> * **Total Variants:** 1,852,266[cite: 1]
> * **Total Samples:** 2,504 (1,270 Females / 1,233 Males / 1 Ambiguous)[cite: 1]
> * **Founders:** 2,497[cite: 1]

### External GWAS Statistics
* **Path:** `\data\external\gwas`[cite: 1]
* **Description:** Type 2 Diabetes (T2D) summary statistics used purely for ranking SNPs by association strength[cite: 1]. 
* **Retained Columns:** Chromosome, Position, EffectAllele, NonEffectAllele, PvalAssociation, Beta, SE, EAF, Ncases, Ncontrols, Neff[cite: 1].

### Downstream Processed Matrices
* **Generation Script:** `\scripts\make_top1000_matrix.py`[cite: 1]
* **Logic:** Matched the external GWAS coordinates against the clean master 1,000 Genomes variant list, selecting the top 1,000 SNPs based on the lowest $p$-values[cite: 1].
* **Generated Outputs:**
  * `\data\processed\selected_snps.txt` (SNP ID List)[cite: 1]
  * `\data\processed\selected_snps_metadata.csv` (Metadata)[cite: 1]
  * `\data\processed\snp_matrix_top1000.raw` (Raw dosage matrix: 2,504 samples $\times$ 1,000 SNPs)[cite: 1]