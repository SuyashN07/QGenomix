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

## Repository Files

### `README.md`

This document describes the QGenomix workflow, dataset expectations, machine-learning pipeline, evaluation metrics, and planned comparisons between classical, biological, and quantum approaches.

### `quantum_blueprint.txt`

Defines the interface contracts between the quantum-processing stages of the project. It specifies expected input shapes, value ranges, and output formats for each quantum layer.

## Classical Machine-Learning Workflow

The classical machine-learning pipeline is designed to perform the following steps:

1. Load the SNP feature matrix from `snp_matrix.npy`.
   - Expected shape: `2504 × 1000`
   - Rows represent individuals.
   - Columns represent SNP features.

2. Load phenotype labels from `epistatic_phenotypes.csv`.
   - Expected number of labels: `2504`
   - Labels are binary values: `0` or `1`.

3. Split the dataset into training, validation, and test subsets.
   - Training set: 70%
   - Validation set: 15%
   - Test set: 15%

4. Train baseline classical models using `baselines.py`.

5. Generate predictions on the validation and test sets.
   - Hard predictions: binary values, `0` or `1`
   - Probability predictions: continuous values between `0.0` and `1.0`

6. Evaluate model performance using `metrics.py`.
   - AUC-ROC
   - F1 score

7. Calculate a biology-driven polygenic risk score using `prs.py`.
   - This component depends on the availability of beta values.
   - The resulting score will also be evaluated using the common metrics pipeline.

8. Compare the performance of:
   - Classical machine-learning models
   - Polygenic risk scores
   - Future quantum-machine-learning models

## Planned Project Structure

The following files are expected as the project develops:

```text
QGenomix/
├── README.md
├── quantum_blueprint.txt
├── snp_matrix.npy
├── epistatic_phenotypes.csv
├── baselines.py
├── metrics.py
└── prs.py
```

## Data Requirements

The expected input data consists of:

- A NumPy SNP matrix containing genetic features.
- A CSV file containing binary phenotype labels.
- Beta values for calculating polygenic risk scores.

The input features and labels must have matching sample counts.

## Evaluation

Models are evaluated using metrics suitable for binary classification:

- **AUC-ROC** measures how well a model separates the two classes across different classification thresholds.
- **F1 score** measures the balance between precision and recall.

Using the same evaluation metrics across all approaches enables a consistent comparison between classical models, polygenic risk scores, and quantum models.

## Quantum-Model Integration

The quantum portion of QGenomix is specified in `quantum_blueprint.txt`. The planned quantum pipeline will generate quantum-derived features that can be passed to a classical machine-learning layer for final prediction.

## Project Status

The repository currently contains the project workflow documentation and the quantum interface specification. Implementation files and datasets can be added as the project progresses.
