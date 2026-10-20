"""Track 2b: do the EFTCAMB python getters change the state of a results object?

In one process per model: read_ini + get_results, then
  E1 = get_eft_functions(a), S1 = get_scale_evolution(k, a)            (a down to 1e-8: forces ThermoData re-init)
  E2, S2 = the same calls again on the same results object
  Cl2 = results.get_cmb_power_spectra() after the getters, compared with Cl1 taken before them
  H2 = results.hubble_parameter(z) after vs before
Expected: E2 == E1, S2 == S1 (excluding the known never-initialised fields EFTcdotdot, EFTGamma1PP, V_scf),
Cl2 == Cl1, H2 == H1, all bit for bit.
Usage: python t2b_getters.py [tree] [model,model,...]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

SKIP = {"EFTcdotdot", "EFTGamma1PP", "V_scf"}
FAST = {"l_max_scalar": "800", "get_tensor_cls": "F", "feedback_level": "0", "EFTCAMB_write_background": "F"}
MODELS = ["2_PEFT_Omega_const_1", "3_RPH_KT_OmegaDE_1", "4_Des_F_R_LCDM_1", "4_Des_mc5e_CPL_1",
          "5_quint_exponential", "5_Kmouflage_3", "5_Horava_eta_1", "5_hdsk_TG_1",
          "008p5-beyond-horndeski"]


def child(tree: str, ini_path: str, workdir: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    w = Path(workdir) / "w.ini"
    w.write_text("".join(f"{k} = {v}\n" for k, v in FAST.items()) + f"DEFAULT({ini_path})\n")
    cp = camb.read_ini(str(w))
    r = camb.get_results(cp)
    z = np.linspace(0, 5, 21)
    cl1 = r.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"].copy()
    h1 = r.hubble_parameter(z)
    a = np.logspace(-8, 0, 60)
    k = [1e-3, 0.05, 0.5]

    def snap():
        _, e = cp.EFTCAMB.get_eft_functions(r, a)
        _, s = cp.EFTCAMB.get_scale_evolution(r, k, a)
        return {f: np.array(e[f]) for f in e.dtype.names}, {f: np.array(s[f]) for f in s.dtype.names}

    e1, s1 = snap()
    e2, s2 = snap()
    cl2 = r.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"]
    h2 = r.hubble_parameter(z)
    bad_e = sorted(f for f in e1 if f not in SKIP and not np.array_equal(e1[f], e2[f], equal_nan=True))
    bad_s = sorted(f for f in s1 if f not in SKIP and not np.array_equal(s1[f], s2[f], equal_nan=True))
    worst_s = max((float(np.nanmax(np.abs(s1[f] - s2[f])) / (np.nanmax(np.abs(s1[f])) or 1.0)) for f in bad_s),
                  default=0.0)
    print("JSON" + json.dumps({"eft_diff_fields": bad_e, "evol_diff_fields": bad_s[:12], "n_evol_diff": len(bad_s),
                               "worst_evol_rel": worst_s, "cl_same": bool(np.array_equal(cl1, cl2)),
                               "H_same": bool(np.array_equal(h1, h2))}), flush=True)


def one(tree: str, model: str) -> str:
    ini = harness.PARAMS_DIR / f"{model}.ini"
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    if not ini.exists():
        ini = HERE.parents[1] / "inis" / f"{model}.ini"
    with harness.slot(1), tempfile.TemporaryDirectory() as d:
        p = subprocess.run([sys.executable, __file__, "_child", str(harness.tree_path(tree)), str(ini), d], cwd=d,
                           env=env, capture_output=True, text=True, timeout=540, preexec_fn=harness._nice)
    js = [ln[4:] for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        return f"{model:26s} FAILED rc={p.returncode} {p.stderr.strip().splitlines()[-1:] if p.stderr else ''}"
    d = json.loads(js[0])
    ok = not d["eft_diff_fields"] and not d["n_evol_diff"] and d["cl_same"] and d["H_same"]
    return (f"{model:26s} {'PASS' if ok else 'FAIL'} eft_diff={d['eft_diff_fields']} n_evol_diff={d['n_evol_diff']} "
            f"worst_evol_rel={d['worst_evol_rel']:.2e} {d['evol_diff_fields']} cl_same={d['cl_same']} "
            f"H_same={d['H_same']}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["_child"]:
        child(*sys.argv[2:5])
    else:
        tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
        models = sys.argv[2].split(",") if len(sys.argv) > 2 else MODELS
        with ThreadPoolExecutor(int(os.environ.get("JOBS", "2"))) as ex:
            for line in ex.map(lambda m: one(tree, m), models):
                print(line, flush=True)
