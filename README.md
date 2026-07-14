# QGenomix Classical ML:

Workflow for the Classical ML part:

1. Load data
   → snp_matrix.npy (2504 people × 1000 SNPs)
   → epistatic_phenotypes.csv (2504 labels: 0 or 1)

2. Split into train/validation/test sets
   → 70% train, 15% val, 15% test

3. baselines.py trains 4 models on the training set

4. Each model predicts on the validation set
   → hard predictions (0 or 1)
   → probabilities (0.0 to 1.0)

5. metrics.py scores those predictions
   → AUC-ROC, F1

6. prs.py (once unblocked) computes a biology-driven score (will be done when the beta values are available)
   → also evaluated by metrics.py

7. Everything gets compared
   → Classical models vs PRS vs (eventually) quantum model
