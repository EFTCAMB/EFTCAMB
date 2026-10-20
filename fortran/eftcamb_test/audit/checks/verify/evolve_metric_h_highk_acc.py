"""Verifier: high-k P(k) discrepancy of evolve_metric_h=T (massless nu, cubic on) vs accuracy boost."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "evolve_metric_h"
INI = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
ML = {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0",
      "accuracy_boost": "2"}
jobs = [(f"ml_acc2_{h}", ML | {"EFTCAMB_evolve_metric_h": h}) for h in ("F", "T")]
with ThreadPoolExecutor(2) as ex:
    st = dict(ex.map(lambda j: (j[0], harness.run_model("A0", INI, OUT / f"{j[0]}.npz", overrides=j[1])), jobs))
print({n: s["status"] for n, s in st.items()})
ks = [0.1, 1.0, 3.0, 10.0, 30.0, 50.0, 100.0]
for a, b in (("ml_acc2_F", "ml_acc2_T"), ("masslessnu_T", "ml_acc2_T"), ("masslessnu_F", "ml_acc2_F")):
    ra, rb = harness.load(OUT / f"{a}.npz"), harness.load(OUT / f"{b}.npz")
    pa = np.interp(np.log(ks), np.log(ra["pk_k"]), np.log(ra["pk"][0]))
    pb = np.interp(np.log(ks), np.log(rb["pk_k"]), np.log(rb["pk"][0]))
    print(f"{b} vs {a}: dP " + " ".join(f"k{k:g}:{x:+.2e}" for k, x in zip(ks, np.exp(pb - pa) - 1)))
