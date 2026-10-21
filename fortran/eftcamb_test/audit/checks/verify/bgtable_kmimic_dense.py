"""Dense-grid EFT functions of 5_Kmimic_1 at model_background_num_points = 1000, 1001, 10000 (A0): look for
localized spikes / singular features that the 401-point harness grid misses."""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable_kmimic"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0",
        "get_transfer": "F"}
FIELDS = ["EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTOmegaPPP", "EFTc", "EFTLambda", "EFTcdot", "EFTLambdadot",
          "EFTGamma1V", "EFTGamma1P", "adotoa", "Hdot", "Hdotdot"]


def child(n: str) -> None:
    import numpy as np

    sys.path.insert(0, str(harness.TREES["A0"]))
    import camb

    w = OUT / f"dense_N{n}.ini"
    ov = {**FAST, "model_background_num_points": n}
    w.write_text("\n".join([f"DEFAULT({harness.PARAMS_DIR / '5_Kmimic_1.ini'})", f"output_root = {OUT / 'dense'}",
                            "EFTCAMB_write_background = F"] + [f"{k} = {v}" for k, v in ov.items()]) + "\n")
    cp = camb.read_ini(str(w), no_validate=True)
    res = camb.get_background(cp)
    a = np.logspace(-8, 0, 200001)
    _, arr = cp.EFTCAMB.get_eft_functions(res, a)
    np.savez(OUT / f"dense_N{n}.npz", a=a, **{f: np.asarray(arr[f]) for f in FIELDS})


def main() -> None:
    import numpy as np

    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    ns = ["1000", "1001", "10000"]
    for n in ns:
        with harness.slot(1):
            r = subprocess.run([sys.executable, __file__, "_child", n], env=env, capture_output=True, text=True,
                               timeout=600)
        if r.returncode:
            print(n, "FAILED", r.stderr[-2000:])
    d = {n: np.load(OUT / f"dense_N{n}.npz") for n in ns}
    a = d["1000"]["a"]
    rows = {}
    for f in FIELDS:
        line = []
        for n in ns:
            y = d[n][f]
            # spike indicator: max |y - running median| / median |y|
            k = 25
            pad = np.pad(y, k, mode="edge")
            from numpy.lib.stride_tricks import sliding_window_view

            med = np.median(sliding_window_view(pad, 2 * k + 1), axis=1)
            spike = np.abs(y - med) / max(np.median(np.abs(y)), 1e-300)
            i = int(np.argmax(spike))
            line.append((n, float(np.max(np.abs(y))), float(a[np.argmax(np.abs(y))]), float(spike[i]), float(a[i]),
                         bool(np.all(np.isfinite(y)))))
        y0, y1 = d["1000"][f], d["1001"][f]
        rel = float(np.max(np.abs(y1 - y0)) / max(np.max(np.abs(y0)), 1e-300))
        rows[f] = dict(per_N=line, rel_diff_1000_1001=rel)
        print(f"{f:12s} rel|N1001-N1000|={rel:9.2e}  " + "  ".join(
            f"N{n}: max={m:9.2e}@a={am:8.2e} spike={s:8.2e}@a={sa:8.2e} fin={fin}" for n, m, am, s, sa, fin in line))
    (OUT / "dense_summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_child":
        child(sys.argv[2])
    else:
        main()
