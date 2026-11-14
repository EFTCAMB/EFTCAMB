"""Sanity check of the EFTCAMB test spectra: every numeric output file must be finite and not identically zero,
and the temperature power spectrum must be positive.

numdiff accepts two files that both contain nan, or that are both all zeros, so a broken run can match a broken
reference: this check runs on the new results before the comparison with the legacy ones.

usage: python check_spectra.py RESULTS_DIR [MODEL ...]   (exit code 1 if a file fails)
"""

import glob
import os
import sys

import numpy as np

# files whose second column is the temperature power spectrum, l (l+1) C_l / 2 pi
TT_FILES = ("_scalCls.dat", "_lensedCls.dat", "_totCls.dat", "_lensedtotCls.dat", "_lenspotentialCls.dat")


def read_numeric(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            try:
                rows.append([float(x.replace("D", "E")) for x in s.split()])
            except ValueError:
                # text lines (e.g. a header without #): a nan or inf token still fails below
                if any(t in s.lower().split() for t in ("nan", "inf", "-inf", "+inf")):
                    return None
    width = min((len(r) for r in rows), default=0)
    return np.array([r[:width] for r in rows]) if rows else np.zeros((0, 0))


def check_file(path):
    data = read_numeric(path)
    if data is None or (data.size and not np.isfinite(data).all()):
        return "nan or inf"
    if data.size == 0:
        return "no numeric data"
    values = data[:, 1:] if data.shape[1] > 1 else data
    if not np.any(values):
        return "all zeros"
    if path.endswith(TT_FILES) and data.shape[1] > 1 and not (data[:, 1] > 0).all():
        return "TT <= 0"
    return None


def main():
    results = sys.argv[1]
    models = sys.argv[2:]
    files = sorted(glob.glob(os.path.join(results, "*.dat")))
    if models:
        files = [f for f in files if any(os.path.basename(f).startswith(m + "_") for m in models)]
    failed = [(f, r) for f in files if (r := check_file(f))]
    for f, reason in failed:
        print(f"  SANITY FAIL: {os.path.basename(f)}: {reason}")
    print(f"  Checked {len(files)} files, {len(failed)} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
