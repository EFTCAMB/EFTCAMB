"""Perturbation evolution of one test model (track t10b-convergence).

python evol.py MODEL OUT.npz [--tree A0] [--k 1e-3,1e-2] [--amin 1e-6] [--na 2000] [--set key=value ...]
Runs in a subprocess importing camb from the tree; uses cheap CAMB settings (no Cls beyond lmax 200, no lensing,
no transfer, no windows) since only the background/thermal history is needed by get_scale_evolution.
Saves evol_<field>[k, a] for all timestep-cache fields.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

CHEAP = {
    "get_tensor_cls": "F",
    "num_redshiftwindows": "0",
    "do_lensing": "F",
    "l_max_scalar": "300",
    "get_transfer": "F",
    "feedback_level": "1",
}


def worker(cfg: dict) -> None:
    sys.path.insert(0, cfg["tree"])
    import camb

    cp = camb.read_ini(cfg["ini"], no_validate=True)
    for k, v in cfg.get("python_set", {}).items():
        obj = cp
        parts = k.split(".")
        for p in parts[:-1]:
            obj = getattr(obj, p)
        setattr(obj, parts[-1], v)
    res = camb.get_results(cp)
    a = np.logspace(np.log10(cfg["amin"]), np.log10(cfg.get("amax", 1.0)), cfg["na"])
    fields, arr = cp.EFTCAMB.get_scale_evolution(res, cfg["k"], a)
    out = {f"evol_{f}": np.ascontiguousarray(arr[f]) for f in arr.dtype.names}
    out["evol_k"] = np.asarray(cfg["k"])
    out["evol_a"] = a
    out["meta"] = np.array(json.dumps(cfg))
    np.savez_compressed(cfg["out"], **out)


def run(model: str, out: str | Path, tree: str = "A0", k=(1e-3,), amin: float = 1e-6, na: int = 2000,
        overrides: dict | None = None, python_set: dict | None = None, threads: int = 1, cheap: bool = True,
        amax: float = 1.0) -> int:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ini = Path(model) if model.endswith(".ini") else harness.PARAMS_DIR / f"{model}.ini"
    ov = dict(CHEAP) if cheap else {}
    ov.update(overrides or {})
    wrapper = out.with_suffix(".wrapper.ini")
    lines = [f"DEFAULT({ini.resolve()})", f"output_root = {out.with_suffix('')}", "EFTCAMB_write_background = F"]
    lines += [f"{a} = {b}" for a, b in ov.items()]
    wrapper.write_text("\n".join(lines) + "\n")
    cfg = {"tree": str(harness.tree_path(tree)), "ini": str(wrapper), "out": str(out), "k": list(map(float, k)),
           "amin": amin, "amax": amax, "na": na, "python_set": python_set or {}}
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with harness.slot(threads), open(out.with_suffix(".log"), "w") as f:
        p = subprocess.run([sys.executable, __file__, "_worker", json.dumps(cfg)], cwd=out.parent, env=env,
                           stdout=f, stderr=subprocess.STDOUT, timeout=3600, preexec_fn=lambda: os.nice(10))
    return p.returncode


def main() -> None:
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        worker(json.loads(sys.argv[2]))
        return
    p = argparse.ArgumentParser()
    p.add_argument("model")
    p.add_argument("out")
    p.add_argument("--tree", default="A0")
    p.add_argument("--k", default="1e-3")
    p.add_argument("--amin", type=float, default=1e-6)
    p.add_argument("--na", type=int, default=2000)
    p.add_argument("--amax", type=float, default=1.0)
    p.add_argument("--set", nargs="*", default=[])
    p.add_argument("--full", action="store_true", help="do not apply the cheap settings")
    a = p.parse_args()
    rc = run(a.model, a.out, a.tree, [float(x) for x in a.k.split(",")], a.amin, a.na,
             dict(s.split("=", 1) for s in a.set), cheap=not a.full, amax=a.amax)
    print("rc", rc)


if __name__ == "__main__":
    main()
