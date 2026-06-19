import argparse
import csv
from pathlib import Path

import numpy as np


def load_metadata(path: Path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            row["idx"] = idx
            row["p"] = float(row["p"])
            row["score"] = float(row["score"])
            rows.append(row)
    return rows


def write_csv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def select_candidate_pairs(metadata, matrix, n_pairs: int, ld_threshold: float):
    # Rank likely useful pairs by GWAS strength. Skip pairs that are too correlated.
    order = sorted(range(len(metadata)), key=lambda i: metadata[i]["score"], reverse=True)
    max_snps_to_scan = min(250, len(order))
    candidate_indices = order[:max_snps_to_scan]

    rows = []
    seen = set()

    for a_pos, i in enumerate(candidate_indices):
        xi = matrix[:, i]
        xi_std = float(np.std(xi))
        if xi_std == 0:
            continue

        for j in candidate_indices[a_pos + 1:]:
            pair_key = (min(i, j), max(i, j))
            if pair_key in seen:
                continue

            xj = matrix[:, j]
            xj_std = float(np.std(xj))
            if xj_std == 0:
                continue

            r = float(np.corrcoef(xi, xj)[0, 1])
            if not np.isfinite(r):
                continue
            r2 = r * r
            if r2 > ld_threshold:
                continue

            score = float(metadata[i]["score"] * metadata[j]["score"])
            rows.append({
                "pair_rank": len(rows) + 1,
                "snp_idx_a": i,
                "snp_idx_b": j,
                "snp_id_a": metadata[i]["id"],
                "snp_id_b": metadata[j]["id"],
                "chrom_a": metadata[i]["chrom"],
                "pos_a": metadata[i]["pos"],
                "chrom_b": metadata[j]["chrom"],
                "pos_b": metadata[j]["pos"],
                "score_a": metadata[i]["score"],
                "score_b": metadata[j]["score"],
                "interaction_score": score,
                "r2": r2,
            })
            seen.add(pair_key)

            if len(rows) >= n_pairs:
                return rows

    return rows


def simulate_phenotypes(matrix, sample_ids, pair_rows, n_truth_pairs: int, prevalence: float, seed: int):
    rng = np.random.default_rng(seed)
    truth_pairs = pair_rows[:n_truth_pairs]

    betas = rng.choice([-1.0, 1.0], size=n_truth_pairs) * rng.uniform(0.55, 1.10, size=n_truth_pairs)
    liability = np.zeros(matrix.shape[0], dtype=np.float64)

    truth_rows = []
    for k, (pair, beta) in enumerate(zip(truth_pairs, betas), start=1):
        i = int(pair["snp_idx_a"])
        j = int(pair["snp_idx_b"])
        interaction = matrix[:, i].astype(np.float64) * matrix[:, j].astype(np.float64)
        liability += beta * interaction
        truth_rows.append({
            "truth_rank": k,
            "snp_idx_a": i,
            "snp_idx_b": j,
            "snp_id_a": pair["snp_id_a"],
            "snp_id_b": pair["snp_id_b"],
            "beta_interaction": float(beta),
        })

    noise = rng.normal(0.0, 1.0, size=matrix.shape[0])
    liability = liability + noise
    threshold = float(np.quantile(liability, 1.0 - prevalence))
    labels = (liability >= threshold).astype(int)

    phenotype_rows = []
    for sample_id, label, score in zip(sample_ids, labels, liability):
        phenotype_rows.append({
            "sample_id": sample_id,
            "case_control": int(label),
            "liability_score": float(score),
        })

    return phenotype_rows, truth_rows, threshold


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--processed-dir", required=True)
    parser.add_argument("--external-dir", required=True)
    parser.add_argument("--n-pairs", type=int, default=500)
    parser.add_argument("--n-truth-pairs", type=int, default=20)
    parser.add_argument("--ld-threshold", type=float, default=0.2)
    parser.add_argument("--prevalence", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    matrix = np.load(args.matrix)
    metadata = load_metadata(Path(args.metadata))

    with Path(args.samples).open("r", encoding="utf-8") as f:
        sample_ids = [line.strip() for line in f if line.strip()][1:]

    if matrix.shape[0] != len(sample_ids):
        raise SystemExit(f"Matrix rows ({matrix.shape[0]}) do not match sample IDs ({len(sample_ids)})")
    if matrix.shape[1] != len(metadata):
        raise SystemExit(f"Matrix columns ({matrix.shape[1]}) do not match metadata rows ({len(metadata)})")

    pair_rows = select_candidate_pairs(metadata, matrix, args.n_pairs, args.ld_threshold)
    if len(pair_rows) < args.n_pairs:
        raise SystemExit(f"Only generated {len(pair_rows)} candidate pairs; expected {args.n_pairs}")

    phenotype_rows, truth_rows, threshold = simulate_phenotypes(
        matrix=matrix,
        sample_ids=sample_ids,
        pair_rows=pair_rows,
        n_truth_pairs=args.n_truth_pairs,
        prevalence=args.prevalence,
        seed=args.seed,
    )

    processed_dir = Path(args.processed_dir)
    external_dir = Path(args.external_dir)

    write_csv(
        processed_dir / "epistasis_candidates.csv",
        ["pair_rank", "snp_idx_a", "snp_idx_b", "snp_id_a", "snp_id_b", "chrom_a", "pos_a", "chrom_b", "pos_b", "score_a", "score_b", "interaction_score", "r2"],
        pair_rows,
    )
    write_csv(
        external_dir / "epistatic_phenotypes.csv",
        ["sample_id", "case_control", "liability_score"],
        phenotype_rows,
    )
    write_csv(
        external_dir / "ground_truth_pairs.csv",
        ["truth_rank", "snp_idx_a", "snp_idx_b", "snp_id_a", "snp_id_b", "beta_interaction"],
        truth_rows,
    )

    case_count = sum(row["case_control"] for row in phenotype_rows)
    print(f"matrix_shape={matrix.shape}")
    print(f"candidate_pairs={len(pair_rows)}")
    print(f"truth_pairs={len(truth_rows)}")
    print(f"cases={case_count}")
    print(f"controls={len(phenotype_rows) - case_count}")
    print(f"case_rate={case_count / len(phenotype_rows):.4f}")
    print(f"liability_threshold={threshold:.6f}")
    print(f"wrote={processed_dir / 'epistasis_candidates.csv'}")
    print(f"wrote={external_dir / 'epistatic_phenotypes.csv'}")
    print(f"wrote={external_dir / 'ground_truth_pairs.csv'}")


if __name__ == "__main__":
    main()
