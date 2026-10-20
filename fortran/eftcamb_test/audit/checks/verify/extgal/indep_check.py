"""Independent check of extgal-adotoa-stale-grhonu.

Route: harness.run_model (full CAMB run, get_eft_functions on the 401-point a grid) and an analytic test at a=1 built
from the ini parameters only (no brentq, no EFTCAMB parameter cache): with S=2 the tracker equation at a=1 is
  x^4 = Ov + Om x^2  ->  x^2 = (Om + sqrt(Om^2 + 4 Ov))/2,
and with Om = 1 - Ov (all matter incl. massive nu) x = 1.  Dropping the massive nu (Om -> 1 - Ov - Onu) gives x < 1.
Also: tracker solution with S=2 has an analytic closed form at every a, compared with eft_adotoa on the grid.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

SCR = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad"
os.environ.setdefault("AUDIT_ROOT", SCR + "/audit")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import harness  # noqa: E402

OUT = Path(SCR) / "audit/work/verify/extgal"
INI = HERE / "008p7-ext-galileon.ini"
OVR = {"EFT_ghost_stability": "F", "EFT_gradient_stability": "F", "EFT_additional_priors": "F"}


def analytic(a, Om_tot, Ov):
    # S = 2: y = (x/a)^2 solves y^2 - Om_tot y - Ov = 0
    y = 0.5 * (Om_tot + np.sqrt(Om_tot**2 + 4 * Ov))
    return a * np.sqrt(y)


CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
for kv in json.loads(sys.argv[3]).items():
    pass
res = camb.get_background(pars)
a = np.logspace(-8, 0, 401)
names, arr = pars.EFTCAMB.get_eft_functions(res, a)
Hdt = a * res.hubble_parameter(1 / a - 1) / 299792.458
i05 = int(np.argmin(abs(a - 0.5)))
print("JSON" + json.dumps(dict(a=a.tolist(), adotoa=arr["adotoa"].tolist(), grhonu_tot=arr["grhonu_tot"].tolist(),
                               Hdt=Hdt.tolist(), i05=i05)))
"""


def bg_run(extra: dict) -> dict:
    import json, subprocess
    OUT.mkdir(parents=True, exist_ok=True)
    w = OUT / "w.ini"
    lines = [f"DEFAULT({INI})", "feedback_level = 0", "EFTCAMB_feedback_level = 0", "EFTCAMB_write_background = F"]
    lines += [f"{k} = {v}" for k, v in {**OVR, **extra}.items()]
    w.write_text("\n".join(lines) + "\n")
    with harness.slot(1):
        r = subprocess.run([sys.executable, "-c", CHILD, str(harness.TREES["A0"]), str(w), "{}"], capture_output=True,
                           text=True, cwd=OUT, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=600)
    js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
    if not js:
        raise SystemExit(r.stdout[-2000:] + r.stderr[-2000:])
    d = json.loads(js[-1][4:])
    return {k: (np.array(v) if isinstance(v, list) else v) for k, v in d.items()}


def main() -> None:
    for tag, extra in [("nu", {}), ("nonu", {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046"})]:
        d = bg_run(extra)
        a = d["a"]
        ad = d["adotoa"]
        gnu = d["grhonu_tot"]
        p = harness.parse_ini(INI)
        h = float(p.get("hubble", 67.0)) / 100
        H0 = h * 100 / 299792.458
        Ob = float(p["ombh2"]) / h**2
        Oc = float(p["omch2"]) / h**2
        Onu = gnu[-1] / (3 * H0**2)  # a=1: grhonu_tot (per a^2 normalisation, a=1) / 3H0^2
        # radiation: from CAMB background (photons+massless nu): Omega_r a^-4 ~ tiny at a>=0.2; use rho at a=1
        Tcmb = float(p.get("temp_cmb", 2.7255))
        Og = 2.4728e-5 * (Tcmb / 2.7255) ** 4 / h**2
        Nmassless = float(p.get("massless_neutrinos", 3.046)) if tag == "nonu" else float(p.get("massless_neutrinos"))
        Or = Og * 7 / 8 * (4 / 11) ** (4 / 3) * Nmassless
        Ov = 1 - Ob - Oc - Og - Or - Onu
        print(f"  params: h={h} Ob={Ob:.5f} Oc={Oc:.5f} Onu={Onu:.5f} Og={Og:.3e} Or={Or:.3e} Ov={Ov:.5f}")
        sel = a >= 0.2
        Om_with = (Ob + Oc) * a**-3 + (Og + Or) * a**-4 + gnu / (3 * H0**2 * a**2)
        Om_wo = (Ob + Oc) * a**-3 + (Og + Or) * a**-4
        ex_w = analytic(a, Om_with, Ov) * H0
        ex_wo = analytic(a, Om_wo, Ov) * H0
        r_w = ad[sel] / ex_w[sel] - 1
        r_wo = ad[sel] / ex_wo[sel] - 1
        print(f"  a>=0.2: max|code/analytic(with nu)-1| = {np.max(np.abs(r_w)):.3e}; "
              f"max|code/analytic(no nu)-1| = {np.max(np.abs(r_wo)):.3e}")
        print(f"  adotoa(a=1)/H0-1 = {ad[-1] / H0 - 1:+.3e}; analytic(a=1, no nu)/H0-1 = {ex_wo[-1] / H0 - 1:+.3e}")
        i = np.searchsorted(a, 0.2)
        print(f"  grid a[{i - 1}]={a[i - 1]:.5f} code/an-1={ad[i - 1] / ex_w[i - 1] - 1:+.3e}; "
              f"a[{i}]={a[i]:.5f} code/an-1={ad[i] / ex_w[i] - 1:+.3e}")
        print(f"  CAMB H(z) [dtauda] vs analytic(with nu): a=0.5 {d['Hdt'][d['i05']]/ex_w[d['i05']]-1:+.3e}, a=1 {d['Hdt'][-1]/ex_w[-1]-1:+.3e}")


if __name__ == "__main__":
    main()
