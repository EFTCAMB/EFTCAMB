"""Check candidate fixes (track 2): reuse of results.Params / get_BAO, and pickle / deepcopy / repr round trips.

For each model, in its own subprocess on the given tree:
  - read_ini + get_background  -> reference H(z), conformal time and EFT functions
  - results.get_BAO(z)         -> H must equal the reference
  - get_background(results.Params) -> must equal the reference
  - copy.deepcopy(cp), pickle round trip, eval(repr(cp)) -> get_background must equal the reference bit for bit
Usage: python t2_fixcheck.py TREE_PATH [model,model,...]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

MODELS = ["2_PEFT_Omega_const_1", "2_PEFT_gamma3_lin_1", "3_RPH_KT_OmegaDE_1", "4_Des_F_R_LCDM_1",
          "4_Des_mc5e_CPL_1", "5_quint_exponential", "5_Kmouflage_3", "5_ADE_1", "5_Horava_eta_1", "5_scg",
          "5_hdsk_TG_1"]


def child(tree: str, model: str, out: str) -> None:
    sys.path.insert(0, tree)
    import copy
    import pickle

    import numpy as np

    import camb
    from camb import model as cmodel  # noqa: F401  (needed by eval(repr))

    z = np.linspace(0, 3, 31)
    a = np.logspace(-4, 0, 50)

    def fp(cp):
        r = camb.get_background(cp)
        _, arr = cp.EFTCAMB.get_eft_functions(r, a)
        return r, np.concatenate([r.hubble_parameter(z), [r.conformal_time(0.0)], arr["EFTOmegaV"], arr["EFTc"],
                                  arr["EFTLambda"], arr["adotoa"]])

    wrapper = harness._wrapper_ini(harness.PARAMS_DIR / f"{model}.ini", Path(out).parent, f"fix_{model}",
                                   {"feedback_level": "0"})
    cp = camb.read_ini(str(wrapper))
    r, ref = fp(cp)
    res = {"model": model}

    def attempt(name, fn):
        try:
            res[name] = fn()
        except Exception as e:  # noqa: BLE001
            res[name] = f"error:{type(e).__name__}:{e}"
        Path(out).write_text(json.dumps(res))
        print(name, res[name], flush=True)

    attempt("get_BAO_H_equal", lambda: bool(np.allclose(r.get_BAO([0.5, 1.0])[:, 1], r.hubble_parameter([0.5, 1.0]),
                                                         rtol=1e-12, atol=0)))
    attempt("reuse_Params_exact", lambda: bool(np.array_equal(fp(r.Params)[1], ref)))
    attempt("deepcopy_exact", lambda: bool(np.array_equal(fp(copy.deepcopy(cp))[1], ref)))
    attempt("pickle_exact", lambda: bool(np.array_equal(fp(pickle.loads(pickle.dumps(cp)))[1], ref)))
    attempt("repr_eval_exact", lambda: bool(np.array_equal(fp(eval(repr(cp), {"camb": camb}))[1], ref)))


def main() -> None:
    tree = sys.argv[1]
    models = sys.argv[2].split(",") if len(sys.argv) > 2 else MODELS
    outdir = harness.AUDIT_ROOT / "work" / "t2-state" / "fixcheck" / Path(tree).name
    outdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    keys = ["get_BAO_H_equal", "reuse_Params_exact", "deepcopy_exact", "pickle_exact", "repr_eval_exact"]
    for m in models:
        out = outdir / f"{m}.json"
        out.unlink(missing_ok=True)
        p = subprocess.run([sys.executable, __file__, "_child", tree, m, str(out)], env=env, capture_output=True,
                           text=True, cwd=outdir)
        res = json.loads(out.read_text()) if out.exists() else {}
        tail = [ln for ln in (p.stderr or "").splitlines() if "error" in ln.lower() or ln.startswith("At line")][-2:]
        cells = " ".join(f"{k}={str(res.get(k, 'not-run'))[:40]}" for k in keys)
        print(f"{m:24s} rc={p.returncode:4d} {cells} {' | '.join(tail)[:200]}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "_child":
        child(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        main()
