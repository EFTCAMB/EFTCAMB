"""Verifier (ic1-pidot-transient-noisy-coefficients): TG_2 single-mode glitch at k = 0.0104/Mpc.
get_scale_evolution (A0, 5_hdsk_TG_2 default) for k around 0.0104: first output time with pi != 0 (the IC is applied at
max(turn-on, mode start time)), pi'/(aH pi) there, and pi at a = 1e-3 and 1e-2, compared with neighbouring k."""

import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "tgmode"
OUT.mkdir(parents=True, exist_ok=True)
CHILD = r'''
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
res = camb.get_results(cp)
ks = [0.00941, 0.00989, 0.0102, 0.0104, 0.0106, 0.01093, 0.01149]
a = np.unique(np.concatenate([np.logspace(-8, -4, 4001), [1e-3, 1e-2, 1e-1]]))
_, arr = cp.EFTCAMB.get_scale_evolution(res, ks, a)
out = {}
for i, k in enumerate(ks):
    pi, pid, H = (np.asarray(arr[f][i], float) for f in ("pi", "pidot", "adotoa"))
    j = int(np.argmax(pi != 0))
    out[str(k)] = {"a_first": a[j], "pi_first": pi[j], "pidot_over_aHpi_first": pid[j] / H[j] / pi[j],
                   "pi_1e-3": float(pi[a == 1e-3][0]), "pi_1e-2": float(pi[a == 1e-2][0]),
                   "pi_1e-1": float(pi[a == 1e-1][0])}
json.dump(out, open(sys.argv[3], "w"), indent=1)
'''
wrapper = harness._wrapper_ini(harness.PARAMS_DIR / "5_hdsk_TG_2.ini", OUT, "tg2_mode", {})
for tag, tree in (("A0", "A0"), ("B1", "B1")):
    res = OUT / f"tg2_mode_{tag}.json"
    if not res.exists():
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1), open(OUT / f"tg2_mode_{tag}.log", "w") as f:
            rc = subprocess.run([sys.executable, "-c", CHILD, str(harness.tree_path(tree)), str(wrapper), str(res)],
                                cwd=OUT, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=500,
                                preexec_fn=harness._nice).returncode
        print(tag, "rc", rc)
    if res.exists():
        print(f"== {tag}")
        for k, v in json.loads(res.read_text()).items():
            print(f"  k={float(k):.5f} IC at a={v['a_first']:.3e}: pi={v['pi_first']:+.3e} pi'/(aH pi)={v['pidot_over_aHpi_first']:+.3e}"
                  f"  pi(1e-3)={v['pi_1e-3']:+.4e} pi(1e-2)={v['pi_1e-2']:+.4e} pi(0.1)={v['pi_1e-1']:+.4e}")
