"""Check 13: small API consistency checks with EFT on (run in one subprocess on tree A0).

(a) quintessence + EFTCAMB_effective_w0wa = T: EFTCAMB.get_effective_w0wa() on the user params (before / after
    get_results) and on results.Params (the initialised copy CAMB used for Halofit).
(b) get_zre_from_tau(params, tau_ini) == Reion.redshift of the run that used tau_ini (GR, hdsk_TG_1, Kmouflage).

usage: python c13_small.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

from t8c_lib import PY, WORK, harness, model_ini

CODE = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
out = {}
cp = camb.read_ini(sys.argv[2])
out["w0wa_user_before"] = cp.EFTCAMB.get_effective_w0wa()
r = camb.get_results(cp)
out["w0wa_user_after"] = cp.EFTCAMB.get_effective_w0wa()
out["w0wa_results_Params"] = r.Params.EFTCAMB.get_effective_w0wa()
out["eft_flag_effective_w0wa"] = cp.EFTCAMB.read_parameters().get("EFTCAMB_effective_w0wa")
for ini in sys.argv[3:]:
    p = camb.read_ini(ini)
    res = camb.get_background(p, no_thermo=False)
    tau = p.Reion.optical_depth
    out[ini.split("/")[-1]] = {"tau": tau, "zre_run": res.Params.Reion.redshift, "zre_from_tau": camb.get_zre_from_tau(p, tau)}
print("JSON" + json.dumps(out))
"""


def main():
    out = WORK / "c13"
    out.mkdir(parents=True, exist_ok=True)
    ov = {"get_tensor_cls": "F", "get_scalar_cls": "F", "get_transfer": "F", "feedback_level": "0",
          "num_redshiftwindows": "0"}
    q = harness._wrapper_ini(model_ini("5_quint_exponential"), out, "quintW", ov | {"EFTCAMB_effective_w0wa": "T"})
    zs = [harness._wrapper_ini(model_ini(m), out, m, ov) for m in ("1_EFT_GR", "5_hdsk_TG_1", "5_Kmouflage_1")]
    env = dict(os.environ, OMP_NUM_THREADS="2")
    p = subprocess.run([PY, "-c", CODE, str(harness.tree_path("A0")), str(q), *map(str, zs)], capture_output=True,
                       text=True, cwd=out, env=env)
    lines = [x for x in p.stdout.splitlines() if x.startswith("JSON")]
    if not lines:
        print(p.stdout[-2000:], p.stderr[-2000:])
        sys.exit(1)
    res = json.loads(lines[-1][4:])
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
