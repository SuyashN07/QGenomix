
Action: file_editor create /app/QGenomix_alternatives_detail.md --file-text "# QGenomix — Alternative Datasets & Pipeline Architectures

**Document Type:** Project Companion Specification
**Scope:** Dataset alternatives, candidate pipeline architectures (Quantum / Classical-ML / Hybrid / Graph / Foundation-Model / Federated / AutoML), reproducible workflow design.
**Status:** Independent specification — companion to the QGenomix Project Plan & Dataset Design document.
**Convention:** Tools are referenced **by name only**; implementation details are deliberately omitted.

---

## 1. Purpose & Positioning

QGenomix targets **polygenic disease risk prediction with quantum-classical hybrid learning**, with explicit attention to **epistatic (gene–gene) interactions** that classical additive models under-capture. This document complements the primary plan by:

1. Curating a **comprehensive catalog of alternative datasets** that can substitute, augment, or stress-test the default 1000 Genomes Phase 3 + simulated phenotype pairing.
2. Enumerating **all credible pipeline topologies** — from pure-classical baselines to quantum-only, hybrid stacks, graph-based, federated, AutoML/NAS, and genomic foundation-model pipelines.
3. Defining a **canonical reproducible workflow** with pluggable backends so experiments across datasets × pipelines × hyperparameters can be launched by a single configuration switch.
4. Eliminating manual orchestration overhead through declarative configs, deterministic seeds, content-addressed artifacts, and CI-validated runs.

This document is intentionally **architecture-and-data-centric**; it does not enumerate team assignments, person-level responsibilities, or sprint-level scheduling.

---

## 2. Objectives

1. Provide **drop-in dataset alternatives** so QGenomix can be benchmarked across population, disease-specific, synthetic, summary-statistic, and multi-omics regimes without code rewrites.
2. Provide a **menu of pipeline topologies** that share a common interface (data → features → model → evaluation → report) so they are mutually substitutable.
3. Make every experiment **reproducible by construction**: identical seeds, pinned environments, content-addressed data, versioned configs, and machine-checked results.
4. Make experimentation **cheap to launch and impossible to lose**: a single command should reproduce any past run, and a single config edit should spawn a new one.
5. Define **clear acceptance criteria** for what qualifies as a valid dataset choice and a valid pipeline implementation.

---

## 3. Scope

**In scope**
- Public and semi-public genomic datasets suitable for polygenic risk / epistasis research.
- Synthetic and simulated datasets to inject known ground-truth interactions.
- Pipeline architectures spanning quantum, classical, hybrid, graph, federated, AutoML, and foundation-model paradigms.
- Reproducibility stack (orchestrators, experiment trackers, environment lockers, data versioners).
- Cross-cutting evaluation protocol (metrics, statistical tests, ablations).

**Out of scope**
- Personnel allocation, task division, or sprint planning.
- Clinical deployment, regulatory submission, or patient-facing UX.
- Detailed tool internals — only the tool **name** is recorded; usage specifics live in code/configs.

---

## 4. Alternative Datasets — Comprehensive Catalog

Datasets are grouped by access tier and intended role. Every entry lists: **role, scale, label availability, epistasis suitability, access tier, and substitution notes**.

### 4.1 Population-Reference Genotype Resources

| # | Dataset | Role | Scale (approx.) | Labels | Epistasis Use | Access | Substitution Notes |
|---|---------|------|------------------|--------|----------------|--------|---------------------|
| 1 | **1000 Genomes Phase 3** | Default reference | 2,504 samples × ~88M variants | None (simulate) | Excellent for simulated epistasis | Open | Current baseline. |
| 2 | **Human Genome Diversity Project (HGDP)** | Population diversity | 929 samples × ~75M variants | None (simulate) | Good — diverse ancestries | Open | Use to stress-test ancestry robustness. |
| 3 | **gnomAD (variant frequency)** | Allele-frequency prior | ~800K samples (aggregate) | None | Indirect — prior for variant filtering | Open | Use for QC/filtering, not as genotype source. |
| 4 | **TOPMed (Bravo summary)** | Frequency reference | ~150K samples (aggregate) | None | Indirect | Open (summary) | Replace gnomAD where deeper coverage needed. |
| 5 | **Simons Genome Diversity Project (SGDP)** | Deep coverage diversity | 300 samples × WGS | None (simulate) | Good — high-quality WGS | Open | Smaller but very clean; useful for unit/integration tests. |

### 4.2 Disease-Cohort & Biobank Resources (Controlled Access)

| # | Dataset | Role | Scale (approx.) | Labels | Epistasis Use | Access | Substitution Notes |
|---|---------|------|------------------|--------|----------------|--------|---------------------|
| 6 | **UK Biobank** | Primary real-world target | ~500K participants × genotype + phenotype | Disease, quantitative traits | Excellent — large N enables interaction power | Application required | The most direct real-data substitute for the simulated phenotype workflow. |
| 7 | **FinnGen** | Founder-population validation | ~500K participants | Disease endpoints | Good — enriched rare variants | Application required | Replication cohort. |
| 8 | **All of Us Research Program** | Diverse U.S. cohort | ~250K+ with genotypes | EHR-linked phenotypes | Good — ancestry diversity | Application (Workbench) | Use for ancestry-stratified evaluation. |
| 9 | **Estonian Biobank** | European replication | ~200K participants | EHR phenotypes | Good | Application | Cross-cohort replication. |
| 10 | **BioBank Japan** | East-Asian ancestry | ~200K participants | 47 disease endpoints | Good | Application | Trans-ancestry transfer tests. |
| 11 | **Million Veteran Program (MVP)** | Large diverse cohort | ~1M participants | EHR phenotypes | Excellent — diversity + scale | Application (VA) | High-statistical-power replication. |

### 4.3 GWAS Summary-Statistic Resources (for PRS & Beta Sourcing)

| # | Dataset | Role | Use | Access |
|---|---------|------|-----|--------|
| 12 | **GWAS Catalog (EBI/NHGRI)** | Curated GWAS index | Source independent betas for PRS | Open |
| 13 | **PGS Catalog** | Pre-validated PRS models | Drop-in PRS baselines, weight files | Open |
| 14 | **OpenGWAS / IEU GWAS DB** | Harmonized summary stats | Beta sourcing, MR-style sanity checks | Open |
| 15 | **FinnGen Summary Stats** | Founder-population betas | Replication PRS | Open |
| 16 | **PAN-UKB Summary Stats** | Multi-ancestry UKB GWAS | Trans-ancestry PRS | Open |

### 4.4 Synthetic & Simulated Datasets (Ground-Truth Epistasis)

| # | Dataset / Simulator Output | Role | Why It Matters |
|---|----------------------------|------|----------------|
| 17 | **HAPGEN2-style simulated cohorts** | Synthetic genotypes preserving LD | Arbitrary sample size; controlled MAF & LD. |
| 18 | **msprime / stdpopsim simulated genomes** | Coalescent-based whole genomes | Demographically realistic; ideal for ablations. |
| 19 | **GAMETES / EpiSim / Triskelion outputs** | Pure-epistasis ground truth | Known interacting SNP pairs/triplets — gold standard for EDR@K. |
| 20 | **Custom additive + interaction injection on 1000G** | Default QGenomix simulator | Current baseline; lives in the dataset-design doc. |
| 21 | **Negative-control permuted phenotypes** | Null distribution | Calibrates false-positive rates of interaction detection. |

### 4.5 Multi-Omics & Auxiliary Modalities (Future Extensions)

| # | Dataset | Modality | Role |
|---|---------|----------|------|
| 22 | **GTEx** | Tissue-specific expression (eQTL) | Functional priors on SNP relevance. |
| 23 | **ENCODE / Roadmap Epigenomics** | Regulatory annotations | Variant prioritization, feature masking. |
| 24 | **ClinVar** | Clinical variant interpretation | Sanity-check top-ranked SNPs. |
| 25 | **Open Targets Genetics** | Variant-to-gene-to-trait | Post-hoc biological interpretation. |
| 26 | **TCGA (germline subset)** | Cancer cohorts | Disease-specific PRS extension. |

### 4.6 Dataset Selection Decision Matrix

| Research Question | Recommended Primary | Recommended Replication | Synthetic Companion |
|-------------------|---------------------|--------------------------|---------------------|
| Method validation, EDR@K benchmarking | 1000 Genomes + injected epistasis | HGDP + injected epistasis | GAMETES |
| Real-disease generalization | UK Biobank | FinnGen / MVP | — |
| Trans-ancestry robustness | All of Us | BBJ + PAN-UKB | HGDP-simulated |
| PRS-only baseline | PGS Catalog weights on UKB | FinnGen | — |
| Rare-variant epistasis | TOPMed-frequency-filtered cohort | FinnGen | stdpopsim rare-variant scenarios |
| Quantum-kernel scaling study | Down-sampled 1000G slices | Synthetic msprime | GAMETES |

### 4.7 Acceptance Criteria for a Dataset to Enter the Pipeline

A dataset qualifies as a usable QGenomix input only if it satisfies **all** of:

1. **Schema-conformant**: convertible to the canonical `(samples × SNP-dosage matrix, phenotype vector, covariates)` triple.
2. **QC-passable**: survives standard variant/sample QC (call-rate, MAF, HWE, heterozygosity, relatedness).
3. **Reproducibly retrievable**: pinned version or checksum; ephemeral URLs are rejected.
4. **License-clear**: open or with documented access approval; controlled-access subsets are isolated behind a separate config.
5. **Ground-truth-or-simulatable phenotype**: either real labels, or a documented simulator with a fixed seed.

---

## 5. Pipeline Architectures — Full Catalog

All pipelines share the **canonical 5-stage contract**:

```
[Stage A] Data Ingest & QC
        ↓
[Stage B] Feature Construction
        ↓
[Stage C] Model Training (pipeline-specific core)
        ↓
[Stage D] Inference & Calibration
        ↓
[Stage E] Evaluation, Statistical Testing & Reporting
```

Each pipeline below is a **drop-in replacement for Stage C** (and sometimes B), preserving Stages A, D, E so results are directly comparable.

### 5.1 Pipeline P0 — Pure Classical-ML Baseline

**Intent:** establish a strong, fast, well-understood reference.

```
SNPs ──► QC ──► Feature Matrix (dosage + PRS + covariates) ──► [LogReg | RF | XGBoost | SVM] ──► Predictions ──► Metrics
```

- **Pros:** fast, mature tooling, strong baselines, interpretable.
- **Cons:** weak on high-order epistasis.
- **Primary role:** lower bound; every other pipeline must beat this.

### 5.2 Pipeline P1 — Pure Quantum (VQC-Only) Classifier

**Intent:** test whether a quantum model alone can predict risk.

```
SNP Pairs ──► Angle Encoding ──► VQC Ansatz ──► Measurement ──► Quantum Logits ──► Predictions
```

- **Pros:** isolates quantum signal; clean ablation reference.
- **Cons:** limited by qubit count and noise; weak when interactions are sparse vs additive bulk.
- **Role:** ablation; quantifies the quantum-only contribution.

### 5.3 Pipeline P2 — Quantum-Kernel-Only (QSVM)

**Intent:** use the quantum circuit purely as a feature map; classify with a classical kernel SVM.

```
SNP Pairs ──► Quantum Feature Map ──► Kernel Matrix K(x,x') ──► Classical SVM ──► Predictions
```

- **Pros:** no variational training; deterministic given the feature map; theoretically grounded.
- **Cons:** O(N²) kernel cost; PCA/landmark approximations needed.
- **Role:** non-variational quantum baseline.

### 5.4 Pipeline P3 — Hybrid Stacking Ensemble (QGenomix Default)

**Intent:** combine quantum features, quantum kernel, and classical learners.

```
                ┌──► VQC Features (N × 1000) ──┐
SNPs ──► QC ──► │                              ├──► Concatenate ──► [QSVM | XGBoost] ──► LogReg Meta-Learner ──► Predictions
                └──► PRS + Covariates ─────────┘
                              │
                              └──► SHAP on XGBoost ──► Adaptive SNP-Pair Re-Selection ──► (loop back to VQC)
```

- **Pros:** strongest expected performance; explainable via SHAP; adaptive loop closes the feedback gap.
- **Cons:** highest implementation complexity; quantum components are the bottleneck.
- **Role:** primary research artifact.

### 5.5 Pipeline P4 — Graph-Based (GNN over LD / Interaction Graph)

**Intent:** treat the genome as a graph (nodes = SNPs, edges = LD or biological interaction) and learn with a GNN.

```
SNPs + LD Graph ──► Graph Construction ──► GNN (GCN / GAT / GraphSAGE) ──► Pooling ──► Predictions
```

- **Pros:** natively captures structured interactions; scales to whole-genome via sparse graphs.
- **Cons:** graph construction is itself a modeling choice; LD ≠ epistasis.
- **Role:** strong non-quantum interaction-aware baseline; also a candidate for **quantum-graph hybrid** (replace one GNN layer with a VQC).

### 5.6 Pipeline P5 — Foundation-Model / Genomic-LLM Pipeline

**Intent:** use pretrained genomic sequence/variant foundation models to produce embeddings, then classify.

```
Variants / Sequence Windows ──► Pretrained Genomic Foundation Model ──► Per-sample Embedding ──► Classical Head (or VQC head) ──► Predictions
```

- **Pros:** leverages massive pretraining; strong on small downstream cohorts; embeddings are reusable.
- **Cons:** model size, license constraints, sequence-window framing for SNPs is non-trivial.
- **Role:** modern strong baseline; the **VQC-head variant** becomes a fourth hybrid.

### 5.7 Pipeline P6 — AutoML / NAS-Driven Pipeline

**Intent:** automate model selection, hyperparameters, and (optionally) quantum-circuit architecture search.

```
Feature Matrix ──► Search Space (classical models, hyperparams, VQC ansätze) ──► AutoML/NAS Controller ──► Best Pipeline ──► Predictions
```

- **Pros:** removes manual tuning; can discover non-obvious circuit ansätze (Quantum Architecture Search).
- **Cons:** compute-hungry; risk of overfitting the search to the validation split.
- **Role:** meta-pipeline that **outputs another pipeline** as its artifact.

### 5.8 Pipeline P7 — Federated / Privacy-Preserving Pipeline

**Intent:** train across institutions without centralizing genotypes.

```
Site 1 Data ─┐
Site 2 Data ─┼─► Local Training (any of P0–P5) ──► Encrypted Gradients / Weights ──► Federated Aggregator ──► Global Model
Site N Data ─┘
```

- **Pros:** unlocks controlled-access biobanks that cannot share raw data; aligns with regulatory reality.
- **Cons:** convergence is slower; quantum components add a coordination layer; needs secure aggregation.
- **Role:** deployment-realistic pipeline for cross-biobank studies.

### 5.9 Pipeline P8 — Hybrid Variant: Quantum-Augmented GNN

**Intent:** combine the graph view (P4) with quantum feature extraction (P3).

```
SNP Graph ──► GNN Encoder ──► Node Embeddings ──► VQC on top-K Nodes ──► Pooling ──► Classifier ──► Predictions
```

- **Role:** novel research direction; treats VQC as a per-cluster interaction probe.

### 5.10 Pipeline P9 — Tensor-Network Surrogate Pipeline

**Intent:** approximate the VQC with a classical tensor-network (MPS/PEPS) representation to test whether the quantum advantage survives classical simulation.

```
SNP Pairs ──► Tensor-Network Feature Map ──► Contraction ──► Features ──► Classical Head ──► Predictions
```

- **Role:** quantum-advantage falsification test; if P9 ≈ P3, the quantum component is not contributing irreducible value.

### 5.11 Pipeline Comparison Matrix

| ID | Pipeline | Quantum Use | Interaction Capture | Compute Cost | Interpretability | Primary Role |
|----|----------|-------------|----------------------|---------------|------------------|---------------|
| P0 | Classical-ML | None | Weak | Low | High | Lower-bound baseline |
| P1 | VQC-only | Full | Medium | Medium | Low | Quantum ablation |
| P2 | Quantum-Kernel-only | Feature map | Medium | Medium-High | Medium | Non-variational quantum baseline |
| P3 | Hybrid Stacking (default) | Features + Kernel | Strong | High | Medium (via SHAP) | Primary artifact |
| P4 | Graph (GNN) | None | Strong | Medium | Medium | Non-quantum interaction baseline |
| P5 | Genomic Foundation Model | Optional head | Medium-Strong | High | Low-Medium | Modern strong baseline |
| P6 | AutoML / NAS | Optional in search | Variable | Very High | Variable | Meta-search |
| P7 | Federated | Any of P0–P5 | Inherits | High (coordination) | Inherits | Privacy-preserving deployment |
| P8 | Quantum-Augmented GNN | Embedded | Very Strong (claimed) | High | Medium | Novel hybrid |
| P9 | Tensor-Network Surrogate | Simulated | Medium | Medium | Medium | Quantum-advantage falsifier |

### 5.12 Pipeline Acceptance Criteria

A pipeline qualifies for inclusion in the experiment grid only if it:

1. Implements the **canonical 5-stage contract** with identical Stage A and Stage E.
2. Consumes the **canonical dataset triple** without bespoke preprocessing.
3. Emits the **canonical results schema** (predictions, probabilities, metrics, calibration, runtime, memory).
4. Runs end-to-end from a **single config file** plus a **single command**.
5. Is **deterministic given a seed** (or documents the source and magnitude of non-determinism).

---

## 6. The Canonical Reproducible Workflow

The goal is: **one config → one command → one fully reproduced run**, regardless of which dataset (Section 4) or pipeline (Section 5) is selected.

### 6.1 Workflow Topology

```
                ┌──────────────────────────────────────────────┐
                │                Experiment Config             │
                │  (dataset, pipeline, seed, hyperparams, env) │
                └────────────────────┬─────────────────────────┘
                                     │
                                     ▼
        ┌────────────────────────────────────────────────────┐
        │  Stage A — Ingest & QC  (dataset-plugin selected)  │
        └────────────────────────┬───────────────────────────┘
                                 ▼
        ┌────────────────────────────────────────────────────┐
        │  Stage B — Feature Construction                    │
        │  (PRS, covariates, SNP-pair selection, embeddings) │
        └────────────────────────┬───────────────────────────┘
                                 ▼
        ┌────────────────────────────────────────────────────┐
        │  Stage C — Pipeline Core (P0…P9 plugin selected)   │
        └────────────────────────┬───────────────────────────┘
                                 ▼
        ┌────────────────────────────────────────────────────┐
        │  Stage D — Inference & Calibration                 │
        └────────────────────────┬───────────────────────────┘
                                 ▼
        ┌────────────────────────────────────────────────────┐
        │  Stage E — Evaluation, Statistical Tests, Report   │
        └────────────────────────┬───────────────────────────┘
                                 ▼
              Artifacts: metrics.json, predictions.parquet,
              calibration.json, figures/, logs/, run-manifest.json
```

### 6.2 Plugin Contracts

Each dataset and each pipeline is implemented as a **plugin** behind a stable interface:

- **Dataset plugin** must return: `(genotype_matrix, phenotype_vector, covariates, metadata)` with a documented schema and a fixed checksum.
- **Pipeline plugin** must implement: `fit(train) → model`, `predict(model, X) → (y_hat, proba)`, `describe() → manifest`.

This contract is what makes the **dataset × pipeline grid** a single command-line sweep.

### 6.3 Configuration Layering

Configurations are layered (most-general → most-specific), and the final effective config is **hash-stamped** into the run manifest:

```
defaults.yaml
   └── dataset.<name>.yaml          (1000g, ukb, hgdp, gametes, …)
         └── pipeline.<id>.yaml     (P0, P1, …, P9)
               └── experiment.yaml  (seed, splits, sweep axes)
                     └── CLI overrides
```

### 6.4 Reproducibility Stack (Names Only)

The recommended **full stack** — adopt as much as the project warrants:

- **Workflow orchestration:** Snakemake, Nextflow
- **Experiment tracking & registry:** MLflow, Weights & Biases
- **Config management:** Hydra, OmegaConf
- **Data & artifact versioning:** DVC, Git-LFS
- **Environment isolation:** Docker, Singularity/Apptainer, Conda/Mamba
- **CI/CD & gated reproducibility:** GitHub Actions
- **Cluster execution:** SLURM
- **Determinism aids:** fixed seeds, pinned library versions, content-addressed dataset checksums

**Orchestrator recommendation:** standardize on **Snakemake** for genomics-style DAGs because it composes naturally with PLINK-style step graphs; keep **Nextflow** as a supported alternative for HPC sites that prefer it. Track every run in **MLflow**; version data and large artifacts with **DVC**; manage configs with **Hydra**; containerize with **Docker** (dev) and **Singularity** (HPC); gate merges with **GitHub Actions**.

### 6.5 Determinism & Provenance Rules

1. **Seed everything** — Python, NumPy, framework RNGs, dataloader shuffles, quantum sampling.
2. **Pin everything** — library versions, container digests, dataset checksums.
3. **Record everything** — the run manifest captures git SHA, container digest, config hash, dataset checksums, hardware, wall-clock, and final metric hashes.
4. **Re-run-on-demand** — given only `run_manifest.json`, the pipeline must reproduce identical artifacts (within documented quantum-sampling tolerance).
5. **Reject silent drift** — CI compares a small reference run's metric hash on every merge.

### 6.6 Experiment Grid (Default Sweep)

The minimum publishable sweep:

```
datasets   = {1000G+inj, HGDP+inj, GAMETES, UKB-subset(if available), msprime}
pipelines  = {P0, P1, P2, P3, P4, P5, P9}             # P6/P7/P8 are optional extensions
seeds      = {0, 1, 2, 3, 4}
splits     = stratified 5-fold CV
```

Total: `|datasets| × |pipelines| × |seeds|` runs, each producing a uniform results record.

### 6.7 Reporting Layer

A single post-hoc step ingests all `metrics.json` files and produces:

- A **leaderboard table** (mean ± CI per dataset × pipeline).
- **Pairwise statistical tests** (DeLong for AUC, McNemar for predictions, bootstrap CIs).
- **Calibration plots** and **decision-curve analyses**.
- **EDR@K** for epistasis-aware pipelines (P1, P2, P3, P8, P9).
- **Ablation deltas** isolating the quantum contribution (P3 − P0, P3 − P4, P3 − P9).
- An auto-generated **Markdown + LaTeX report** committed alongside the run.

---

## 7. Evaluation Protocol (Shared Across All Pipelines)

- **Classification metrics:** AUC-ROC, AUC-PR, F1, balanced accuracy, Brier score, ECE (calibration).
- **Epistasis-aware metrics:** EDR@K, top-K precision/recall against known interacting pairs (only meaningful on synthetic ground-truth datasets).
- **Statistical tests:** DeLong (AUC), McNemar (predictions), bootstrap 95% CI on metric differences.
- **Robustness checks:** ancestry-stratified metrics, label-permutation null, downsampling curves.
- **Compute metrics:** wall-clock per stage, peak memory, qubit-count, circuit-depth, kernel-matrix build time.

---

## 8. Risk Register (Architecture-Level)

| # | Risk | Likelihood | Impact | Mitigation |
|---|------|-------------|--------|-------------|
| 1 | Controlled-access cohorts unavailable in time | Medium | High | Default to 1000G + simulated phenotypes; keep biobank pipelines behind a feature flag. |
| 2 | Quantum kernel O(N²) blow-up at scale | High | High | PCA reduction, landmark/Nyström approximation, sub-sampling protocols. |
| 3 | Non-determinism from quantum sampling | Medium | Medium | Document tolerance bands; fix shot counts; report sampling-noise CIs. |
| 4 | Foundation-model license / weight access | Medium | Medium | Keep P5 optional; substitute open alternatives if needed. |
| 5 | AutoML overfitting the search split | Medium | High | Nested CV; report search-budget sensitivity. |
| 6 | Federated coordination overhead | High | Medium | Keep P7 as a follow-on milestone, not core. |
| 7 | Reproducibility decay (lib upgrades) | High | High | Lockfiles + container digests + CI reference run. |
| 8 | Quantum-advantage illusion (P3 ≈ P9 or P4) | Medium | Critical | Always run P9 (tensor-network surrogate) and P4 (GNN) as falsifiers. |

---

## 9. Deliverables of This Specification

1. **Dataset Catalog** (Section 4) — 26 entries across population, biobank, summary-stat, synthetic, and multi-omics tiers, with a decision matrix and acceptance criteria.
2. **Pipeline Catalog** (Section 5) — 10 pipeline topologies (P0–P9) with a shared 5-stage contract and comparison matrix.
3. **Reproducible Workflow Blueprint** (Section 6) — plugin contracts, layered configs, determinism rules, default sweep, reporting layer.
4. **Evaluation Protocol** (Section 7) — uniform metrics + statistical tests applicable to every pipeline.
5. **Risk Register** (Section 8) — architecture-level risks and mitigations.

---

## 10. Acceptance Criteria for This Document

This specification is considered complete when:

- [ ] Every alternative dataset has a documented role, scale, access tier, and substitution note.
- [ ] Every pipeline (P0–P9) is described with its data flow, pros/cons, and role.
- [ ] The canonical workflow is fully described and reproducible via a single config + command pattern.
- [ ] The reproducibility stack is named (not detailed) and an orchestrator is recommended.
- [ ] An evaluation protocol applicable across all pipelines is defined.
- [ ] The document contains **no task assignments** and **no tool internals**.

---

## 11. Glossary (Brief)

- **SNP** — Single-Nucleotide Polymorphism.
- **PRS** — Polygenic Risk Score.
- **Epistasis** — non-additive interaction between two or more genetic variants.
- **VQC** — Variational Quantum Circuit.
- **Quantum Kernel** — kernel function evaluated via overlap of quantum-encoded states.
- **EDR@K** — Epistasis Detection Rate among the top-K ranked variant pairs.
- **GNN** — Graph Neural Network.
- **NAS** — Neural Architecture Search (here extended to Quantum Architecture Search).
- **Federated Learning** — distributed training without raw-data exchange.
- **Tensor-Network Surrogate** — classical approximation of a quantum state (e.g., MPS/PEPS) used to test whether quantum-only gains survive classical simulation.

---

*End of specification. Companion document to the QGenomix Project Plan and the QGenomix Dataset Design document.*
"
Observation: Create successful: /app/QGenomix_alternatives_detail.md