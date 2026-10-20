"""Track 2a helper: per-field exact comparison of two thread runs (t1 vs t8 or repeats), listing every
eft_* / evol_* field that is not bit identical, with max relative diff and a magnitude summary of the values."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

D = harness.AUDIT_ROOT / "work" / "t2-state" / "threads"


def diff_fields(pa: Path, pb: Path) -> dict[str, tuple[float, float, float]]:
    A, B = harness.load(pa), harness.load(pb)
    out = {}
    for k in sorted(A):
        if not (k.startswith("eft_") or k.startswith("evol_")) or k not in B:
            continue
        x, y = np.asarray(A[k]), np.asarray(B[k])
        if x.dtype.kind not in "fc" or x.shape != y.shape:
            continue
        if np.array_equal(x, y, equal_nan=True):
            continue
        scale = float(np.nanmax(np.abs(y))) if np.isfinite(y).any() else 0.0
        d = float(np.nanmax(np.abs(x - y)))
        out[k] = (d / scale if scale else float("inf"), float(np.nanmax(np.abs(x))), scale)
    return out


if __name__ == "__main__":
    tree, a, b = sys.argv[1], sys.argv[2], sys.argv[3]
    models = sys.argv[4].split(",") if len(sys.argv) > 4 else None
    d = D / tree
    names = sorted({p.name[: -len(f"_{a}.npz")] for p in d.glob(f"*_{a}.npz")})
    for m in names:
        if models and m not in models:
            continue
        pa, pb = d / f"{m}_{a}.npz", d / f"{m}_{b}.npz"
        if not pb.exists():
            continue
        res = diff_fields(pa, pb)
        txt = ", ".join(f"{k}:{r:.1e}(|a|max={ma:.1e},|b|max={mb:.1e})" for k, (r, ma, mb) in res.items())
        print(f"{m:28s} {txt or 'ALL eft/evol fields bit identical'}")
