QGenomix Mem1 Data Handoff Package
Date: 27 May 2026

Contents:
- snp_matrix.npy: standardized SNP matrix, shape 2504 x 1000
- sample_ids.csv: sample IDs matching rows of snp_matrix.npy
- selected_snps.txt: selected top 1000 SNP IDs
- selected_snps_metadata.csv: SNP metadata and GWAS scores
- epistasis_candidates.csv: 500 SNP-pair candidates for VQC
- epistatic_phenotypes.csv: synthetic case/control labels
- ground_truth_pairs.csv: 20 true simulated SNP-pair interactions
- make_top1000_matrix.py: script used to create SNP matrix
- simulate_phenotypes_and_pairs.py: script used to create labels and pairs

Important:
These are synthetic T2D-like labels, not real diabetes diagnoses.
Rows in snp_matrix.npy match sample_ids.csv and epistatic_phenotypes.csv.
