"""Adversarial check for eft-par-cache-realloc-abort: is "same result on reuse" actually the right expectation?

For each model: fresh full calculation from the ini parameters, then a second full calculation
camb.get_results(results.Params.copy()) on a new CAMBdata (State%CP carries the RGR-overwritten
EFTCAMB_pert_turn_on, the already-initialised model object and the filled eft_par_cache).
Compares lensed TT/EE/TE (TE / sqrt(TT EE)), sigma8 and EFTCAMB_pert_turn_on between the two runs.
On A0 the second calculation is expected to abort (the finding); on the fixed tree it should match.

Usage: python verify_par_cache_reuse_spectra.py <tree> [model ...]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
OUT = ROOT / "work" / "verify" / "par_cache_realloc" / "reuse_spectra"


def child(tree: str, ini: str, out: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    pars = camb.read_ini(ini)
    t_user = float(pars.EFTCAMB.EFTCAMB_pert_turn_on) if hasattr(pars, "EFTCAMB") else None
    r1 = camb.get_results(pars)
    cl1 = r1.get_total_cls(2500)
    s81 = float(r1.get_sigma8_0())
    t_state = float(r1.Params.EFTCAMB.EFTCAMB_pert_turn_on)
    res = {"turn_on_user": t_user, "turn_on_state_after_first": t_state, "sigma8_1": s81}
    Path(out).write_text(json.dumps(res))
    np.save(out + ".cl1.npy", cl1)
    print("FIRST_DONE", flush=True)
    p2 = r1.Params.copy()
    r2 = camb.get_results(p2)
    cl2 = r2.get_total_cls(2500)
    res.update(sigma8_2=float(r2.get_sigma8_0()), turn_on_state_after_second=float(r2.Params.EFTCAMB.EFTCAMB_pert_turn_on))
    np.save(out + ".cl2.npy", cl2)
    Path(out).write_text(json.dumps(res))
    print("SECOND_DONE", flush=True)


def compare(out: Path) -> dict:
    import numpy as np

    res = json.loads(out.read_text())
    c2 = Path(str(out) + ".cl2.npy")
    if not c2.exists():
        return res
    a, b = np.load(str(out) + ".cl1.npy"), np.load(c2)
    ell = slice(2, 2501)
    d = {}
    for i, name in ((0, "TT"), (1, "EE")):
        d[name] = float(np.max(np.abs(b[ell, i] / a[ell, i] - 1)))
    d["TE"] = float(np.max(np.abs(b[ell, 3] - a[ell, 3]) / np.sqrt(a[ell, 0] * a[ell, 1])))
    res["max_rel_cl"] = d
    res["sigma8_rel"] = res["sigma8_2"] / res["sigma8_1"] - 1
    res["bit_identical_cl"] = bool(np.array_equal(a, b))
    return res


def main() -> None:
    tree = sys.argv[1]
    tree = tree if os.path.isabs(tree) else str(ROOT / "trees" / tree)
    models = sys.argv[2:] or ["4_Des_F_R_LCDM_1", "2_PEFT_Omega_lin_1", "5_Kmouflage_1"]
    label = Path(tree).name
    OUT.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    summary = {}
    for m in models:
        with tempfile.NamedTemporaryFile("w", suffix=".ini", dir=OUT, delete=False) as f:
            f.write(f"DEFAULT({PARAMS / (m + '.ini')})\noutput_root = {OUT}/{label}_{m}\n"
                    "EFTCAMB_write_background = F\nfeedback_level = 0\nEFTCAMB_feedback_level = 0\n")
            ini = f.name
        out = OUT / f"{label}_{m}.json"
        for p in (out, Path(str(out) + ".cl1.npy"), Path(str(out) + ".cl2.npy")):
            p.unlink(missing_ok=True)
        with harness.slot(1):
            p = subprocess.run(["nice", "-n", "10", sys.executable, __file__, "_child", tree, ini, str(out)],
                               capture_output=True, text=True, env=env, cwd=str(OUT), timeout=900)
        err = [ln for ln in (p.stdout + p.stderr).splitlines() if "runtime error" in ln or ln.startswith("At line")]
        r = compare(out) if out.exists() else {}
        r.update(rc=p.returncode, err=[e.replace(str(ROOT), "$S") for e in err])
        summary[m] = r
        print(f"== tree={label} model={m}: {json.dumps(r)}", flush=True)
    (OUT / f"summary_{label}.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_child":
        child(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        main()
