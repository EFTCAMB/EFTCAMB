"""Repro: EFTCAMB returns (and NaN-checks) never-written memory in three timestep-cache fields.

1) V_scf of Horndeski models with Horndeski_designer_V = F (e.g. 5_hdsk_TG_1, 5_hdsk_TG_3):
   fortran/eftcamb/08f_full_models/008p0_Horndeski.f90:1778  eft_cache%V_scf = self%Vscf%value(...)
   but self%Vscf%y is only filled when designerV is true (008p0_Horndeski.f90:2996-3000); the table is allocated
   by EquispacedLinearIntepolateFunction1DInit and never written, so get_eft_functions returns heap garbage.
2) EFTcdotdot and EFTGamma1PP: EFTCAMBTimestepCacheInit (03_abstract_EFTCAMB_cache.f90:402-586) resets every
   field except these two; they are only set by compute_additional_derivs (stability check, mass stability only),
   so get_scale_evolution returns uninitialised stack memory for them (all full-mapping models, designer f(R)).
3) The same fields are tested by EFTCAMBTimestepCacheIsNan inside EFTTestStability (09_EFTCAMB_stability.f90:207)
   before they are ever written (first sampled point has a < EFTCAMB_pert_turn_on), so the stability verdict
   depends on uninitialised memory (see the -finit-derived/NaN-table instrumented run quoted in the finding).

V_scf: the same model is computed 4 times in one process (get_results and get_background alternated); EFTcdotdot/EFTGamma1PP: two separate
processes. Expected: identical, deterministic values.
Run: /Users/marco/Software/miniforge3/bin/python uninit-cache-fields-garbage.py [tree]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
PARAMS = Path(__file__).resolve().parents[2] / "parameters"
CHEAP = {"get_transfer": "F", "get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100",
         "feedback_level": "0", "EFTCAMB_feedback_level": "0", "EFTCAMB_write_background": "F"}


def child(tree: str, model: str, what: str, workdir: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    ini = Path(workdir) / f"{model}_wrap.ini"
    ini.write_text("\n".join(f"{k} = {v}" for k, v in CHEAP.items()) + f"\nDEFAULT({PARAMS / (model + '.ini')})\n")
    cp = camb.read_ini(str(ini))
    a = np.array([1e-3, 1e-2, 0.1, 0.5, 1.0])
    out = {}
    if what == "eft_functions":
        # the same model four times in one process: the never-written V_scf table picks up whatever the
        # heap held (typically numbers left by an earlier calculation)
        out["V_scf"], out["phi_scf"] = [], []
        for i in range(4):
            cp = camb.read_ini(str(ini))
            r = camb.get_results(cp) if i % 2 == 0 else camb.get_background(cp)
            _, arr = cp.EFTCAMB.get_eft_functions(r, a)
            out["V_scf"].append([float(x) for x in arr["V_scf"]])
            out["phi_scf"].append([float(x) for x in arr["phi_scf"]])
    else:
        r = camb.get_results(cp)
        _, arr = cp.EFTCAMB.get_scale_evolution(r, [0.01, 0.1], a)
        for f in ("EFTcdotdot", "EFTGamma1PP", "EFTcdot"):
            out[f] = [repr(float(x)) for x in np.asarray(arr[f]).ravel()]
    print("JSON" + json.dumps(out), flush=True)


def run(tree: str, model: str, what: str) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory() as d:
        p = subprocess.run([sys.executable, __file__, "_child", tree, model, what, d], env=env, cwd=d,
                           capture_output=True, text=True)
    for ln in p.stdout.splitlines():
        if ln.startswith("JSON"):
            return json.loads(ln[4:])
    raise RuntimeError(f"child failed rc={p.returncode}: {p.stderr[-800:]}")


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else "A0"
    tree = arg if os.path.isabs(arg) else str(ROOT / "trees" / arg)
    print(f"tree: {tree}")
    bad = 0
    for model in ("5_hdsk_TG_1", "5_hdsk_TG_3"):
        r1 = run(tree, model, "eft_functions")
        print(f"{model}: get_eft_functions V_scf at a=[1e-3,1e-2,0.1,0.5,1] (designer_V = F, table never filled),"
              " get_results/get_background alternated 4x in one process")
        for i, v in enumerate(r1["V_scf"]):
            print(f"   run {i + 1}: V_scf = {['%.6g' % x for x in v]}   phi_scf = {['%.6g' % x for x in r1['phi_scf'][i]]}")
        same = all(v == r1["V_scf"][0] for v in r1["V_scf"])
        print(f"   V_scf identical in the 4 runs: {same};  phi_scf identical: "
              f"{all(v == r1['phi_scf'][0] for v in r1['phi_scf'])}")
        bad += not same
    for model in ("5_quint_exponential", "4_Des_F_R_LCDM_1"):
        r1, r2 = run(tree, model, "evolution"), run(tree, model, "evolution")
        print(f"{model}: get_scale_evolution k=[0.01,0.1], 5 times (first 4 values shown)")
        for f in ("EFTcdotdot", "EFTGamma1PP", "EFTcdot"):
            same = r1[f] == r2[f]
            print(f"   {f:12s} p1={r1[f][:4]}  p2={r2[f][:4]}  identical={same}")
            bad += (not same) and f != "EFTcdot"
    print("expected: identical deterministic values in both processes (0 or the physical value);")
    print(f"observed: {bad} field(s) differ between two identical runs (uninitialised memory)" if bad else
          "observed: all identical")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_child":
        child(*sys.argv[2:6])
    else:
        main()
