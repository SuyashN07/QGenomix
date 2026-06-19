import argparse, csv, gzip, io, math, zipfile, struct
from pathlib import Path

def open_gwas(path):
    path = Path(path)
    if path.suffix.lower() == ".zip":
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if not n.endswith("/") and not n.lower().endswith((".pdf", ".docx"))]
        name = max(names, key=lambda n: z.getinfo(n).file_size)
        return io.TextIOWrapper(z.open(name), encoding="utf-8", errors="replace")
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def norm(s):
    return s.lower().replace("#", "").replace("_", "").replace("-", "").replace(".", "")

def find_col(header, names):
    nh = [norm(x) for x in header]
    for name in names:
        if name in nh:
            return nh.index(name)
    return None

def split_line(line, tab_mode):
    return line.rstrip("\n").split("\t") if tab_mode else line.strip().split()

def write_npy_float32(path, rows, shape):
    magic = b"\x93NUMPY"
    version = b"\x01\x00"
    header = "{'descr': '<f4', 'fortran_order': False, 'shape': %s, }" % (shape,)
    header_bytes = header.encode("latin1")
    pad = 16 - ((10 + len(header_bytes) + 1) % 16)
    header_bytes += b" " * pad + b"\n"
    with open(path, "wb") as f:
        f.write(magic + version + struct.pack("<H", len(header_bytes)) + header_bytes)
        for row in rows:
            f.write(struct.pack("<" + "f" * len(row), *row))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pvar", required=True)
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n", type=int, default=1000)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    pos_to_pvar = {}
    with open(args.pvar, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            chrom, pos, vid, ref, alt = parts[:5]
            chrom = chrom.replace("chr", "")
            pos_to_pvar.setdefault((chrom, pos), []).append((vid, chrom, pos, ref, alt))

    print(f"Loaded {sum(len(v) for v in pos_to_pvar.values())} variants from merged PVAR")

    selected = {}
    with open_gwas(args.gwas) as f:
        header_line = ""
        for line in f:
            if line.strip():
                header_line = line
                break

        tab_mode = "\t" in header_line
        header = split_line(header_line, tab_mode)

        chr_i = find_col(header, ["chr", "chrom", "chromosome", "chromsome"])
        pos_i = find_col(header, ["pos", "bp", "position", "basepairlocation"])
        p_i = find_col(header, ["pvalassociation", "pvalueassociation", "p", "pvalue", "pval", "pvaluefrequentist", "pvaluegc"])

        if chr_i is None or pos_i is None or p_i is None:
            raise SystemExit(f"Could not detect chr/pos/p columns. Header was:\n{header}")

        for line in f:
            if not line.strip():
                continue
            parts = split_line(line, tab_mode)
            if len(parts) <= max(chr_i, pos_i, p_i):
                continue

            chrom = parts[chr_i].replace("chr", "")
            pos = parts[pos_i]
            try:
                p = float(parts[p_i])
            except ValueError:
                continue

            key = (chrom, pos)
            if key not in pos_to_pvar:
                continue

            if p <= 0:
                p = 1e-300
            score = -math.log10(p)

            for vid, c, bp, ref, alt in pos_to_pvar[key]:
                old = selected.get(vid)
                if old is None or p < old["p"]:
                    selected[vid] = {"id": vid, "chrom": c, "pos": bp, "ref": ref, "alt": alt, "p": p, "score": score}

    ranked = sorted(selected.values(), key=lambda r: r["p"])[: args.n]
    if len(ranked) < args.n:
        raise SystemExit(f"Only matched {len(ranked)} SNPs. Stop and inspect GWAS format/build.")

    with open(out / "selected_snps.txt", "w", encoding="utf-8") as f:
        for r in ranked:
            f.write(r["id"] + "\n")

    with open(out / "selected_snps_metadata.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "chrom", "pos", "ref", "alt", "p", "score"])
        w.writeheader()
        w.writerows(ranked)

    print(f"Selected top {len(ranked)} SNPs")
    print("Now converting PLINK .raw to standardized .npy")

    with open(args.raw, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()
        meta_cols = 6
        variant_cols = header[meta_cols:]

        samples = []
        matrix = []
        for line in f:
            if not line.strip():
                continue
            row = line.strip().split()
            if len(row) < meta_cols + 1:
                continue
            samples.append(row[1])
            vals = []
            for x in row[meta_cols:]:
                vals.append(None if x in ("NA", ".", "") else float(x))
            matrix.append(vals)

    n_rows = len(matrix)
    n_cols = len(matrix[0])

    cols = list(zip(*matrix))
    medians, means, stds = [], [], []

    for col in cols:
        present = sorted(x for x in col if x is not None)
        mid = len(present) // 2
        med = present[mid] if len(present) % 2 else (present[mid - 1] + present[mid]) / 2
        filled = [med if x is None else x for x in col]
        mean = sum(filled) / len(filled)
        var = sum((x - mean) ** 2 for x in filled) / len(filled)
        std = math.sqrt(var) or 1.0
        medians.append(med)
        means.append(mean)
        stds.append(std)

    standardized = []
    for row in matrix:
        standardized.append([
            ((medians[j] if row[j] is None else row[j]) - means[j]) / stds[j]
            for j in range(n_cols)
        ])

    write_npy_float32(out / "snp_matrix.npy", standardized, (n_rows, n_cols))

    with open(out / "sample_ids.csv", "w", encoding="utf-8") as f:
        f.write("sample_id\n")
        for s in samples:
            f.write(s + "\n")

    print(f"Wrote {out / 'snp_matrix.npy'} with shape ({n_rows}, {n_cols})")
    print(f"Wrote {out / 'selected_snps_metadata.csv'}")
    print(f"Wrote {out / 'sample_ids.csv'}")

if __name__ == "__main__":
    main()


