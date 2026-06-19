# QGenomix — Full Pipeline Architecture
**Project:** Quantum-Enhanced Epistasis Detection for T2D-like Complex Traits  
**Dataset Mem1:** 1000 Genomes Phase 3 · 2504 samples · Top-1000 GWAS SNPs · 500 epistasis candidates · 20 ground-truth pairs  
**Date:** June 2026  

---

## Part I — Dataset Construction Verification

### 1. Source Data & Provenance

| Layer | Item | Status |
|---|---|---|
| Raw genotypes | 1000 Genomes Phase 3, chr1–22, PLINK2 pgen format (downloaded via `run_download_qc_prune.ps1`) | Confirmed |
| GWAS summary stats | T2D GWAS from IEU OpenGWAS / DIAGRAM consortium | Used for p-value scoring |
| Sample roster | `phase3_corrected.psam` → 2504 unrelated individuals | Confirmed |

### 2. Processing Chain Verification

**Step 1 — Per-chromosome QC (PLINK2)**  
Filters applied: MAF > 0.01, HWE p > 1×10⁻⁶, genotype missingness < 5%, sample missingness < 10%, biallelic SNPs only (`--max-alleles 2 --snps-only just-acgt`), duplicate ID removal.

**Step 2 — LD Pruning (PLINK2)**  
Window=50 SNPs, step=5, r²<0.2 (`--indep-pairwise 50 5 0.2`). Pruned files written per-chromosome to `data/interim/qc_by_chr/`.

**Step 3 — Merge across chromosomes**  
All 22 pruned chromosome files merged into a single PLINK2 dataset at `data/interim/merged/`.

**Step 4 — GWAS-guided top-1000 SNP selection (`make_top1000_matrix.py`)**  
- Reads merged `.pvar` to build chr:pos → variant-id map  
- Reads T2D GWAS summary file; extracts chr, pos, p-value columns with flexible header detection  
- Ranks by –log₁₀(p), takes top-1000 SNPs  
- Converts PLINK `.raw` dosage file to standardized float32 `.npy` matrix: missing values imputed with column median, then z-score standardised per SNP  
- Outputs: `selected_snps.txt`, `selected_snps_metadata.csv`, `snp_matrix.npy`, `sample_ids.csv`

**Step 5 — Candidate pair selection (`simulate_phenotypes_and_pairs.py`)**  
- Scans top-250 highest-scoring SNPs pairwise  
- Rejects pairs with r² > 0.2 (low-LD filter)  
- Scores pairs by product of individual –log₁₀(p) scores  
- Outputs top-500 candidate pairs to `epistasis_candidates.csv`

**Step 6 — Synthetic phenotype simulation**  
- Selects 20 truth pairs (rank 1–20 from candidates)  
- Assigns each a random multiplicative interaction beta ∈ [±0.55, ±1.10]  
- Builds additive epistatic liability: `L = Σ βₖ · xᵢ · xⱼ + ε`, where ε ~ N(0,1)  
- Thresholds at 90th percentile (10% case prevalence)  
- Outputs: `epistatic_phenotypes.csv` (case/control + liability), `ground_truth_pairs.csv` (20 true interacting pairs)

### 3. Dataset Integrity Check

| Check | Expected | Verified |
|---|---|---|
| SNP matrix shape | (2504, 1000) | 2504 rows × 1000 cols, float32 |
| Sample alignment | `snp_matrix.npy` rows ↔ `sample_ids.csv` | Same 2504 samples |
| Phenotype alignment | `epistatic_phenotypes.csv` rows ↔ `sample_ids.csv` | Same 2504 samples |
| Candidate pairs | 500 pairs from top-250 SNPs, r²<0.2 | 500 rows in `epistasis_candidates.csv` |
| Ground truth pairs | 20 pairs, all from rank 1–20 candidates | 20 rows in `ground_truth_pairs.csv` |
| Top SNP score | SNP `6:20679310:T:G`, p=2.43×10⁻³²¹ → –log₁₀p≈320.6 | Matches metadata row 1 |
| LD filter on truth pairs | All 20 pairs have r² < 0.2 | Verified from `epistasis_candidates.csv` |
| Interaction structure | All 20 truth pairs share `snp_idx_a=0` (top GWAS SNP as hub) | Star-topology epistasis |

> **Key structural note:** The synthetic phenotype is built with a single hub SNP (`6:20679310:T:G`, idx=0) interacting with 20 partners. This star-topology is realistic for regulatory-region SNPs but should be noted as a potential model bias when evaluating detector performance — real epistasis is often more dispersed.

### 4. Known Limitations of Current Dataset

- **Synthetic phenotype only** — not real T2D diagnoses. Downstream models are validated against a known ground truth, making them benchmark-compatible but not clinically transferable without real-disease replication.  
- **Star topology** — 20 truth pairs all involve SNP idx=0. This makes recall metrics partially dependent on whether models can identify the hub SNP.  
- **European-ancestry dominance** — 1KG Phase 3 is globally diverse but ~500 of 2504 samples are European. Fine-mapped GWAS hits will have highest power in this group.  
- **n=2504** — borderline for many quantum circuit approaches that are already sample-limited; QSVM and VQC are well-suited to this range.

---

## Part II — Alternative / Supplemental Datasets

The following datasets are structurally compatible with the QGenomix pipeline and can substitute or augment the current 1KG-based dataset.

### Tier 1 — Direct Drop-In Replacements (same genotype + phenotype format)

| Dataset | Samples | Ancestry | Phenotype | Notes |
|---|---|---|---|---|
| **UK Biobank (UKBB)** | ~500,000 | Primarily European | Real T2D ICD-10 codes | Gold standard for T2D GWAS; requires access application |
| **All of Us Research Program** | >650,000 | Diverse US cohort | Real EHR phenotypes | NIH-funded; broad ancestry; T2D well-represented |
| **FinnGen** | ~520,000 | Finnish isolate | Real T2D diagnoses | Very high statistical power in isolate; public summary stats |
| **HUNT Cohort** | ~56,000 | Norwegian | Real T2D | Well-phenotyped; public GWAS results available |
| **DIAGRAM Consortium** | Meta-analysis ~180K cases | Multi-ancestry | Real T2D GWAS | Summary stats used in current pipeline; can extend to individual-level |

### Tier 2 — Genotype-Only (phenotype must be simulated or linked)

| Dataset | Samples | Ancestry | Availability |
|---|---|---|---|
| **TOPMed** | >180,000 | Diverse (8 populations) | dbGaP application required |
| **gnomAD v4** | ~807,000 exomes / ~76,000 genomes | Multi-ancestry | Aggregated allele freq only; not individual-level |
| **HGDP + 1KG (combined)** | ~4,000 | 54 populations globally | Freely available; extends ancestry representation |
| **Simons Diversity Project** | 279 | 142 populations | Ultra-diverse; freely available |
| **African Genome Variation Project** | 320 | 18 African populations | Adds underrepresented ancestry crucial for trans-ethnic epistasis |

### Tier 3 — Benchmark Epistasis Datasets (purpose-built for interaction testing)

| Dataset | Format | Use |
|---|---|---|
| **PLINK epistasis benchmark sets (Urbanowicz et al.)** | Synthetic PLINK | Standard ML epistasis benchmark; known interaction structures |
| **GAMETES-generated datasets** | Synthetic | Controllable epistasis complexity, prevalence, heritability |
| **ATHENA challenge datasets** | Synthetic BED | Used in ATHENA neural network epistasis study; publicly archived |
| **BioBank Japan (BBJ) T2D GWAS** | Summary stats | Non-European validation; ~191K samples |

### Recommendation

For the next iteration of QGenomix, the recommended path is:
1. **Validate** current pipeline against synthetic benchmark sets (GAMETES/Urbanowicz) for ground-truth calibration
2. **Extend** to FinnGen or HUNT for real-disease phenotype under similar ancestry constraints  
3. **Scale** to UKBB for production-level sample size once quantum circuit approaches are stable  
4. **Diversify** by integrating HGDP+1KG for ancestry-aware analysis

---

## Part III — Full Pipeline Architecture

### Overview

The QGenomix pipeline is organized around **five strategy families**, all sharing the same data ingestion and evaluation backbone, making it trivially interchangeable to swap strategies. A master experiment runner wraps everything.

```
DATA INGESTION
    │
    ▼
FEATURE ENGINEERING
    │
    ├──► STRATEGY A: Classical Statistical
    ├──► STRATEGY B: Classical ML / Deep Learning  
    ├──► STRATEGY C: Quantum-Native (VQC / QSVM / QAOA)
    ├──► STRATEGY D: Hybrid Quantum-Classical
    └──► STRATEGY E: Generative / Explainability Models
          │
          ▼
    EVALUATION & COMPARISON
          │
          ▼
    RANKED PAIR OUTPUT → REPORT
```

---

### Stage 0 — Data Ingestion & Validation

**Inputs:**
- `snp_matrix.npy` — (2504 × 1000) float32, z-score standardized
- `sample_ids.csv` — 2504 sample identifiers
- `selected_snps_metadata.csv` — SNP positions, GWAS p-values, –log₁₀p scores
- `epistasis_candidates.csv` — 500 pre-selected low-LD SNP pairs with pair scores
- `epistatic_phenotypes.csv` — case/control labels (0/1) and liability scores
- `ground_truth_pairs.csv` — 20 known true interacting pairs (for evaluation only)

**Operations:**
- Schema validation (shape, dtype, row alignment)
- Train/test split stratified on `case_control` label (e.g., 80/20 or 5-fold CV)
- Optionally: ancestry-aware split using sample metadata
- Cache to HDF5 or parquet for fast repeated access

**Dataflow:**
```
[snp_matrix.npy]               →  DataLoader
[epistatic_phenotypes.csv]     →  LabelLoader
[epistasis_candidates.csv]     →  CandidatePairLoader
                   ↓
           UnifiedDataBundle
           (X_train, X_test, y_train, y_test, candidate_pairs, ground_truth)
```

---

### Stage 1 — Feature Engineering

All strategies receive a **FeatureBundle** derived from the UnifiedDataBundle. The feature engineering stage is modular and switchable via config.

#### 1A — Pair-Level Feature Vectors
For each of the 500 candidate pairs (i,j):
- `x_i`, `x_j` — raw SNP dosage vectors (2504-dim each)
- `x_ij = x_i ⊗ x_j` — elementwise product (interaction feature, 2504-dim)
- `[x_i, x_j, x_ij]` — concatenated 3×2504 matrix = one sample per pair approach
- GWAS scores: `score_i`, `score_j`, `score_i × score_j`, `r²`

#### 1B — Sample-Level Feature Matrices (per-sample, all 1000 SNPs)
- Raw z-score matrix X (2504×1000)
- Pairwise interaction matrix for candidate pairs only: X_pairs (2504×500) where col k = x_i · x_j for pair k
- PCA-reduced: top-k principal components of X (for quantum encoding compression)
- Kernel matrix: pre-computed RBF or quantum kernel (2504×2504) for kernel methods

#### 1C — Graph Features
- LD-based adjacency graph of 1000 SNPs (edge if r² > threshold)
- GWAS-score node features
- Pair-score edge features on candidate pair subgraph

#### 1D — Quantum Encoding Preprocessors
- Amplitude encoding: normalize each SNP pair feature vector to unit norm
- Angle encoding: map features to rotation angles [0, 2π]
- ZZ-feature map kernel precomputation
- IQP (Instantaneous Quantum Polynomial) feature map preparation

---

### Strategy A — Classical Statistical Epistasis Detection

**Goal:** Establish strong classical baselines using established bioinformatics methods.

**Dataflow:**
```
CandidatePairs (500)
    │
    ├──[Logistic Regression]──► β_i, β_j, β_ij interaction test
    │                           Chi-square or LRT p-value per pair
    │
    ├──[PLINK2 --epistasis]───► Exhaustive chi-square interaction scan
    │                           (can extend beyond 500 pairs)
    │
    ├──[BOOST (Bayesian)]────► Bayesian epistasis measure per pair
    │
    └──[REGENIE/SAIGE]───────► Mixed-model interaction test
                                (accounts for population structure)
    │
    ▼
Ranked pair list by interaction p-value
    │
    ▼
Evaluation: AUROC, precision@k, recall@20 against ground truth
```

**Models:**
- Logistic regression with interaction term: `logit(y) = α + β_i·xᵢ + β_j·xⱼ + β_ij·xᵢ·xⱼ`
- ANOVA-based F-test for epistasis
- MDR (Multifactor Dimensionality Reduction)
- BOOST epistasis score

**Hyperparameters:** Interaction term regularization (L1/L2), p-value threshold for reporting.

---

### Strategy B — Classical ML / Deep Learning

**Goal:** Learn non-linear interaction patterns from genotype-phenotype data.

#### B1 — Tree-Based Models

**Dataflow:**
```
X_pairs (2504 × 500 interaction features) + y
    │
    ├──[Random Forest]──────► Feature importance → pair ranking
    ├──[XGBoost/LightGBM]──► SHAP interaction values → pair ranking
    └──[ExtraTrees]────────► Out-of-bag importance ranking
    │
    ▼
Top-k pairs by importance score → Evaluation
```

#### B2 — Neural Epistasis Models

**Dataflow:**
```
Per-pair feature vectors [x_i, x_j, x_ij] (dim=3N)
    │
    ├──[MLP Epistasis Detector]
    │    Input: [x_i; x_j; x_ij] → FC(512)→ReLU→FC(256)→ReLU→FC(1)→sigmoid
    │    Loss: BCE(y_pred, y_case_control)
    │    Gradient saliency → pair interaction score
    │
    ├──[Siamese Network]
    │    Two branches encode x_i and x_j separately → shared FC → interaction head
    │    Learns symmetric interaction score
    │
    ├──[Transformer Epistasis]
    │    SNPs as tokens; self-attention over top-250 SNPs per sample
    │    Attention weights → pair interaction scores
    │    Architecture: 4-head, 2-layer transformer; embed_dim=64
    │
    └──[Variational Autoencoder (VAE)]
         Encode X → latent z → decode; latent space interaction analysis
         Pair score = ||dz/dx_i · dz/dx_j||
    │
    ▼
Pair ranking by gradient/attention/latent score → Evaluation
```

#### B3 — Graph Neural Networks

**Dataflow:**
```
LD graph (1000 nodes, edges where r²>0.05)
Node features: [GWAS score, MAF, chrom_embedding]
Edge features: [r², pair interaction score if in candidate list]
    │
    ├──[GraphSAGE / GCN]──► Node embeddings → pair cosine similarity
    ├──[GAT]──────────────► Attention scores on edges → interaction probability
    └──[GraphTransformer]──► Global attention over all SNP pairs
    │
    ▼
Edge-level interaction scores → Evaluation
```

---

### Strategy C — Quantum-Native Approaches

**Goal:** Exploit quantum superposition, entanglement, and interference to detect non-linear epistatic interactions that may be exponentially expensive classically.

#### C1 — Variational Quantum Classifier (VQC)

**Dataflow:**
```
Per-pair feature vector → Angle/Amplitude Encoding → |ψ(x)>
                                    │
                         Parameterized Ansatz U(θ)
                         (RealAmplitudes, EfficientSU2, TwoLocal)
                                    │
                         Measurement: <Z₀> → binary prediction
                                    │
                         Loss: BCE(pred, y_case_control)
                         Optimizer: COBYLA / SPSA / Adam (parameter-shift)
                                    │
                         Learned θ* → VQC score per pair
    │
    ▼
Top-k pairs by VQC confidence → Evaluation
```

**Circuit variants to explore:**
- `RealAmplitudes` ansatz (shallow, hardware-friendly): 1–3 reps
- `EfficientSU2` ansatz: captures 2-qubit entanglement patterns
- `TwoLocal` with CNOT entanglement: full-qubit connectivity
- `NLocal` with custom gate sets for specific entanglement topologies

**Encoding strategies:**
- ZZFeatureMap: maps pairs to Hilbert space via ZZ interactions (captures multiplicative epistasis naturally)
- PauliFeatureMap with ZXZ terms: richer encoding for 3-way interactions
- Amplitude encoding (log₂(N) qubits for N features): memory-efficient, depth-expensive

**Number of qubits:** 4–10 qubits for pair-level encoding; 10–20 for sample-level PCA-compressed features.

#### C2 — Quantum Support Vector Machine (QSVM / Quantum Kernel Estimation)

**Dataflow:**
```
Pair feature vectors → Quantum Kernel K(xi, xj) = |<ψ(xi)|ψ(xj)>|²
                                    │
                         ZZFeatureMap circuit → kernel matrix (N×N)
                                    │
                         Classical SVM with quantum kernel
                         C parameter tuning via grid search
                                    │
                         Decision function values → pair scores
    │
    ▼
Pair ranking by SVM margin → Evaluation
```

**Key design choices:**
- Feature map reps (1 vs 2): higher reps → richer kernel but more noise
- Kernel alignment vs direct SVM training
- Shot-based estimation vs statevector simulation

#### C3 — Quantum Approximate Optimization Algorithm (QAOA) for Epistasis Pair Selection

**Problem formulation:**  
Frame pair selection as a combinatorial optimization: choose a subset S of k pairs from 500 candidates that maximizes a "epistasis coverage score" — a QUBO (Quadratic Unconstrained Binary Optimization) problem.

**QUBO formulation:**
```
min  -Σᵢ interaction_score(i)·zᵢ  +  λ·Σᵢ≠ⱼ overlap(i,j)·zᵢ·zⱼ
     zᵢ ∈ {0,1}
```
where `overlap(i,j)=1` if pairs share a SNP (penalize redundancy), `λ` is a diversity regularizer.

**Dataflow:**
```
Interaction scores + overlap matrix → QUBO matrix Q
                    │
         QAOA circuit: p alternating layers
         Cost Hamiltonian Hc = Σ Qᵢⱼ ZᵢZⱼ
         Mixer Hamiltonian Hm = Σ Xᵢ
                    │
         Optimize γ,β parameters (COBYLA/SPSA)
                    │
         Sample from |QAOA(γ*,β*)> → binary string z*
                    │
         Selected pairs z* = 1
    │
    ▼
QAOA-selected pair subset → Evaluation
```

**p-layers:** 1–5; higher p → better approximation ratio, more circuit depth.

#### C4 — Quantum Walk-Based Feature Exploration

**Concept:** Model the 500-pair candidate graph as a quantum walk. Pair score nodes become amplitudes. Quantum interference naturally amplifies high-scoring interaction clusters.

**Dataflow:**
```
Candidate pair graph → Adjacency matrix A
                    │
         Quantum walk operator: U = exp(-iAt)
                    │
         Probabilities |<j|U|ψ₀>|² for each pair node
                    │
         High-probability nodes = quantum-amplified pair candidates
    │
    ▼
Walk-amplified ranking → Evaluation
```

#### C5 — Quantum Fisher Information Epistasis Scorer

**Concept:** Quantify how much information SNP pair interactions contribute to phenotype prediction in quantum-encoded space, via the Quantum Fisher Information Matrix (QFIM).

**Dataflow:**
```
Feature-encoded state |ψ(x_pair)>
    │
    QFIM eigenspectrum analysis
    High eigenvalues → high informative subspace
    │
    Trace(QFIM)_pair → pair-level quantum epistasis score
    │
    ▼
QFIM-ranked pairs → Evaluation
```

---

### Strategy D — Hybrid Quantum-Classical

**Goal:** Use quantum circuits as feature transformers inside classical ML pipelines, capturing the best of both worlds.

#### D1 — Quantum Kernel + Classical SVM/Kernel Ridge Regression

**Dataflow:**
```
Features → ZZFeatureMap → Quantum kernel matrix K_Q (N×N)
                │
        K_combined = α·K_Q + (1-α)·K_RBF   (kernel blending)
                │
        Classical SVM / KRR on K_combined
                │
        Pair scores via dual variable weights
    │
    ▼
Hybrid kernel ranking → Evaluation
```

#### D2 — VQC Feature Extractor + Classical Head

**Dataflow:**
```
Pair feature vector → VQC encoder (fixed or trainable)
                            │
                    Intermediate quantum state → measurement bitstring
                            │
                    Bitstring → classical embedding vector
                            │
                    Classical MLP head: FC(128)→ReLU→FC(1)→sigmoid
                            │
                    End-to-end training: quantum+classical parameters
    │
    ▼
Hybrid VQC-MLP pair scores → Evaluation
```

#### D3 — Quantum Attention for Transformer Epistasis

**Concept:** Replace classical dot-product attention in a transformer with a quantum inner-product circuit. SNP tokens attend to each other via quantum overlap.

**Dataflow:**
```
SNP token embeddings (top-250 SNPs) → quantum Q,K,V encoding
                    │
         Attention(Q,K) = |<φ(Q)|φ(K)>|²  (quantum overlap)
                    │
         Quantum-weighted value aggregation
                    │
         Classical FFN + prediction head
    │
    ▼
Quantum-attention pair interaction scores → Evaluation
```

#### D4 — Quantum Boltzmann Machine (QBM) for Epistasis

**Concept:** Model the joint distribution of SNP pairs and phenotype as a Quantum Boltzmann Machine — a transverse-field Ising model where thermal sampling captures correlations.

**Dataflow:**
```
Binary-encoded SNP pairs + case/control label
                    │
         QBM: H = -Σ Jᵢⱼ ZᵢZⱼ - Γ Σ Xᵢ   (transverse-field Ising)
                    │
         Train: maximize log-likelihood via gradient estimation
         (D-Wave or simulated quantum annealing sampler)
                    │
         Learned Jᵢⱼ coupling strengths → epistasis scores
    │
    ▼
QBM coupling-ranked pairs → Evaluation
```

#### D5 — Quantum Generative Adversarial Network (QGAN) for Epistasis Discovery

**Concept:** A QGAN where the quantum generator samples synthetic SNP-pair patterns conditioned on case/control status. Epistasis is detected by comparing generator parameter landscapes between case and control conditions.

**Dataflow:**
```
[Quantum Generator G(z,θ_G)]  ↔  [Classical Discriminator D(x,φ_D)]
        │                                       │
   Latent z + case/ctrl label              Real SNP pair data
        │
   Generator learns p(x_pair | case) and p(x_pair | ctrl)
        │
   KL-divergence between case and ctrl generator landscapes → epistasis signal
    │
    ▼
QGAN-detected interaction patterns → Evaluation
```

#### D6 — Variational Quantum Eigensolver (VQE) for Epistasis Hamiltonian

**Concept:** Encode the epistasis detection problem as finding the ground state of an Ising Hamiltonian where coupling terms Jᵢⱼ represent SNP-pair interaction strengths. VQE finds the minimum-energy configuration that corresponds to the most predictive pair selection.

**Dataflow:**
```
Pair interaction matrix → Ising Hamiltonian H = Σᵢⱼ Jᵢⱼ ZᵢZⱼ
                    │
         VQE: minimize <ψ(θ)|H|ψ(θ)>
         Ansatz: UCCSD-inspired or hardware-efficient
                    │
         Ground state configuration → selected pair subset
    │
    ▼
VQE-selected epistatic pairs → Evaluation
```

---

### Strategy E — Generative, Explainability & Emerging Approaches

#### E1 — Shapley Value Epistasis Attribution

**Dataflow:**
```
Trained model (any strategy) + X_pairs
    │
    SHAP TreeExplainer or KernelSHAP
    Shapley interaction values φᵢⱼ per feature pair
    │
    φᵢⱼ > threshold → epistatic contribution
    │
    ▼
SHAP-attributed pair scores → Evaluation
```

#### E2 — Causal Inference (PC Algorithm / FCI)

**Dataflow:**
```
X_pairs + y
    │
    Constraint-based causal graph learning (PC / FCI)
    Conditional independence tests (Fisher-Z, chi-square)
    │
    Causal graph: directed edges X_i → X_j → Y
    Collider structures: X_i → Y ← X_j  (epistasis marker)
    │
    ▼
Causal graph epistasis edges → Evaluation
```

#### E3 — Diffusion Model for Epistasis Pattern Generation

**Concept:** Train a denoising diffusion probabilistic model (DDPM) on the SNP pair feature space. The score function (diffusion gradient) captures the data distribution geometry and identifies high-density interaction regions.

**Dataflow:**
```
SNP pair features (500 × feature_dim)
    │
    DDPM training: forward noise process + reverse denoising
    │
    Score function ∇log p(x_pair)
    │
    Regions of high score ↔ epistasis-enriched feature space
    │
    ▼
Diffusion score-ranked pairs → Evaluation
```

#### E4 — Topological Data Analysis (TDA) for Epistasis

**Concept:** Compute persistent homology of the SNP pair feature manifold. Topological features (Betti numbers, persistence diagrams) reveal interaction structure not visible in linear statistics.

**Dataflow:**
```
X_pairs embedded in feature space
    │
    Vietoris-Rips filtration → persistence diagram
    │
    Persistent 1-cycles → SNP pairs that form interaction loops
    │
    Persistence lifetime → epistasis interaction strength
    │
    ▼
TDA persistence-ranked pairs → Evaluation
```

#### E5 — Federated Quantum Learning (FQL) — Future Direction

**Concept:** For multi-cohort epistasis discovery without sharing raw genotype data. Each site trains a local quantum model; quantum state transfer encodes the gradient updates, preserving privacy.

**Dataflow:**
```
Site A (1KG)      → Local VQC θ_A
Site B (FinnGen)  → Local VQC θ_B
Site C (BBJ)      → Local VQC θ_C
        │
        Federated averaging: θ_global = mean(θ_A, θ_B, θ_C)
        Quantum secure aggregation protocol
        │
        ▼
Global model → cross-cohort epistasis discovery
```

---

## Part IV — Unified Reproducible Master Pipeline

The master pipeline ties all strategies together under a single configuration-driven runner. No manual workflow is needed — changing a YAML config key switches the entire strategy.

### Directory Structure

```
QGenomix_mem1_step1/
├── data/
│   ├── raw/                     # 1KG Phase 3 pgen files (per-chr)
│   ├── interim/
│   │   ├── qc_by_chr/           # Per-chr QC + pruned pgen
│   │   └── merged/              # Merged pgen across all chromosomes
│   ├── processed/               # snp_matrix.npy, metadata CSVs
│   └── external/                # Phenotypes, ground truth, GWAS files
├── handoff_mem1_data_package/   # Frozen handoff snapshot (read-only)
├── scripts/
│   ├── make_top1000_matrix.py
│   └── simulate_phenotypes_and_pairs.py
├── pipeline/
│   ├── config/
│   │   ├── base.yaml            # Data paths, splits, eval settings
│   │   ├── strategy_A.yaml      # Classical statistical config
│   │   ├── strategy_B.yaml      # Classical ML/DL config
│   │   ├── strategy_C.yaml      # Quantum-native config
│   │   ├── strategy_D.yaml      # Hybrid quantum-classical config
│   │   └── strategy_E.yaml      # Generative/explainability config
│   ├── stages/
│   │   ├── ingest.py            # Stage 0: data loading + validation
│   │   ├── features.py          # Stage 1: feature engineering modules
│   │   ├── train.py             # Stage 2: model training dispatch
│   │   ├── evaluate.py          # Stage 3: evaluation + metrics
│   │   └── compare.py           # Stage 4: multi-strategy comparison
│   ├── models/
│   │   ├── classical/           # Logistic regression, RF, XGBoost, MLP, GNN
│   │   ├── quantum/             # VQC, QSVM, QAOA, QBM, VQE modules
│   │   └── hybrid/              # Hybrid kernel, VQC+MLP, quantum attention
│   ├── run_experiment.py        # Master runner: accepts --config flag
│   └── compare_all.py           # Loads all results, generates comparison report
├── results/
│   └── {strategy}/{run_id}/     # Metrics, predictions, ranked pairs
├── notebooks/
│   └── exploration/             # Interactive analysis notebooks
├── QGenomix_Project_Plan.pdf
├── QGenomix_Pipeline_Architecture.md   (this file)
└── run_download_qc_prune.ps1
```

### Master Config Schema (base.yaml)

```yaml
project: QGenomix
run_id: mem1_v1                     # auto-tagged with timestamp

data:
  snp_matrix: data/processed/snp_matrix.npy
  sample_ids: data/processed/sample_ids.csv
  metadata: data/processed/selected_snps_metadata.csv
  candidates: data/processed/epistasis_candidates.csv
  phenotypes: data/external/epistatic_phenotypes.csv
  ground_truth: data/external/ground_truth_pairs.csv

split:
  method: stratified_kfold          # or train_test, leave_one_out
  n_folds: 5
  test_fraction: 0.20
  random_seed: 42

features:
  mode: pair_level                  # or sample_level, graph, quantum_encoded
  include_gwas_scores: true
  include_r2: true
  pca_components: null              # set int to enable PCA compression

strategy: strategy_C                # A, B, C, D, or E

evaluation:
  metrics: [auroc, auprc, precision_at_k, recall_at_20, f1]
  k_for_precision: [5, 10, 20, 50]
  ground_truth_pairs: 20

output:
  save_ranked_pairs: true
  save_model: true
  save_metrics: true
```

### Experiment Runner (run_experiment.py) — Pseudocode

```python
def run_experiment(config):
    # Stage 0: Ingest
    data = ingest(config.data, config.split)
    
    # Stage 1: Features
    features = build_features(data, config.features)
    
    # Stage 2: Train
    model = dispatch_model(config.strategy, config)
    results = model.fit_predict(features)
    
    # Stage 3: Evaluate
    metrics = evaluate(results.ranked_pairs, config.evaluation)
    
    # Stage 4: Save
    save_results(results, metrics, config.output, config.run_id)
    
    return metrics
```

**Run a single strategy:**
```
python pipeline/run_experiment.py --config pipeline/config/strategy_C.yaml
```

**Run all strategies and compare:**
```
python pipeline/compare_all.py --strategies A B C D E --base-config pipeline/config/base.yaml
```

---

### Stage 3 — Evaluation Framework

All strategies are evaluated identically on the same held-out test set.

#### Primary Metrics

| Metric | Definition | Target |
|---|---|---|
| **AUROC** | Area under ROC curve for pair classification (true vs false epistasis) | > 0.80 = good, > 0.90 = excellent |
| **AUPRC** | Area under precision-recall curve | > 0.50 given 4% base rate (20/500) |
| **Recall@20** | Fraction of 20 ground truth pairs in top-20 ranked predictions | 1.0 = perfect recovery |
| **Precision@k** | Fraction of top-k predictions that are true pairs | k = 5, 10, 20, 50 |
| **Interaction F1** | Harmonic mean of precision and recall at optimal threshold | |

#### Secondary Metrics

| Metric | Purpose |
|---|---|
| **Ranking correlation (Spearman ρ)** | Agreement between predicted and true interaction strength rank |
| **Hub SNP Recovery Rate** | Whether idx=0 (known hub) is in top-10 predicted interactors |
| **Calibration (Brier score)** | Probability calibration quality |
| **Wall-clock runtime** | Practical runtime comparison (quantum sim vs classical) |
| **Circuit depth / gate count** | For quantum strategies: resource estimation |
| **Noise robustness** | Performance degradation under simulated shot noise (1024 vs 8192 shots) |

#### Comparison Dashboard Output

```
Strategy   | AUROC | AUPRC | Recall@20 | Prec@10 | Runtime  | Notes
-----------+-------+-------+-----------+---------+----------+------------------
A-LogReg   | 0.XX  | 0.XX  | XX%       | XX%     | <1s      | Baseline
A-BOOST    | 0.XX  | 0.XX  | XX%       | XX%     | ~10s     |
B-XGBoost  | 0.XX  | 0.XX  | XX%       | XX%     | ~30s     |
B-Transf   | 0.XX  | 0.XX  | XX%       | XX%     | ~5min    |
C-VQC      | 0.XX  | 0.XX  | XX%       | XX%     | ~30min*  | *simulated
C-QSVM     | 0.XX  | 0.XX  | XX%       | XX%     | ~1hr*    | kernel matrix
C-QAOA     | 0.XX  | 0.XX  | XX%       | XX%     | ~20min*  | combinatorial
D-QK+SVM   | 0.XX  | 0.XX  | XX%       | XX%     | ~45min*  | hybrid kernel
D-VQC+MLP  | 0.XX  | 0.XX  | XX%       | XX%     | ~1hr*    | hybrid e2e
E-SHAP     | 0.XX  | 0.XX  | XX%       | XX%     | ~5min    | post-hoc
E-TDA      | 0.XX  | 0.XX  | XX%       | XX%     | ~15min   | topological
```

---

### Stage 4 — Reproducibility Guarantees

1. **Seeded randomness**: all random operations use `seed=42` (configurable); seeds propagated to NumPy, PyTorch, Qiskit Aer simulators.
2. **Immutable handoff snapshot**: `handoff_mem1_data_package/` is treated as read-only; pipeline reads but never writes to it.
3. **Artifact hashing**: input file SHA-256 hashes logged at ingestion; alerts if input data changes unexpectedly.
4. **Environment lock**: `requirements.txt` + `environment.yml` pinned with exact versions.
5. **Run tracking**: each `run_id` is timestamped and logged to `results/`; no results are overwritten.
6. **Config-as-code**: every experiment is fully reproducible by replaying its YAML config.

---

## Part V — Quantum Hardware Considerations

| Aspect | Simulation (current) | Real Hardware (future) |
|---|---|---|
| **Backend** | Qiskit Aer statevector / qasm simulator | IBM Quantum, IonQ, Quantinuum |
| **Qubits** | Up to ~30 (statevector), ~50 (qasm) | 27–127 qubits (IBM Eagle/Heron) |
| **Noise model** | None (ideal) or Aer noise model from device | Real decoherence, gate errors |
| **Shot budget** | Configurable (1024–65536) | 1024–8192 typical per circuit |
| **Circuit depth budget** | Unlimited in simulation | < 100 gates for reliable execution |
| **Recommended first target** | IBM Quantum 27-qubit Falcon (free tier) | |

**Transition strategy:** All quantum circuits should be designed to run within 50-gate depth for hardware compatibility.

---

## Part VI — Experiment Roadmap

```
PHASE 1: Baseline Establishment (Weeks 1-2)
  ├── Run Strategy A (classical statistics) → establish recall@20 floor
  ├── Run Strategy B (XGBoost + MLP) → establish ML ceiling
  └── Validate evaluation framework against ground_truth_pairs.csv

PHASE 2: Quantum Exploration (Weeks 3-6)
  ├── C1: VQC with ZZFeatureMap (4-8 qubits, pair-level)
  ├── C2: QSVM kernel matrix on 500 pairs
  ├── C3: QAOA for pair subset selection (p=1,2,3)
  └── C4: Quantum walk amplification

PHASE 3: Hybrid Integration (Weeks 7-9)
  ├── D1: Quantum kernel blending (alpha sweep)
  ├── D2: VQC encoder + classical MLP head
  └── D5: QBM coupling analysis

PHASE 4: Advanced Models (Weeks 10-12)
  ├── B3: Graph Neural Networks on LD graph
  ├── E2: Causal discovery (FCI algorithm)
  └── E4: Topological Data Analysis

PHASE 5: Comparison & Reporting (Weeks 13-14)
  ├── Run compare_all.py across all strategies
  ├── Noise robustness study for quantum strategies
  ├── Hardware transpilation estimates
  └── Final report + next memory generation handoff
```

---

## Summary

This document establishes:

1. **Dataset verified** — The Mem1 dataset is correctly constructed from 1KG Phase 3 with sound QC, LD pruning, GWAS-guided selection, and synthetic phenotype simulation. The 20 ground-truth epistatic pairs form a star topology around the top GWAS SNP and are suitable for benchmarking all strategies.

2. **Alternative datasets identified** — Ranging from direct drop-in replacements (UKBB, FinnGen) to benchmark sets (GAMETES, ATHENA) and diverse ancestry augmentation (HGDP, AGVP).

3. **Five pipeline strategies** — Classical statistical, Classical ML/DL, Quantum-native, Hybrid quantum-classical, and Generative/Explainability — each with detailed dataflow, model specifications, and hyperparameter guidance.

4. **Unified reproducible master pipeline** — Config-driven experiment runner with standardized evaluation, reproducibility guarantees, and a comparison dashboard for systematic strategy benchmarking.
