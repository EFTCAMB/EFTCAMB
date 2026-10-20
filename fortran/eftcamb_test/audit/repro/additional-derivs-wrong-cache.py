"""Repro: compute_additional_derivs (designer and full-mapping base classes) returns wrong second/third derivatives.

Code: fortran/eftcamb/06p2_abstract_EFTCAMB_designer.f90:61-257 and fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90:169-404.
Inside the finite-difference helper `dummy_helper(temp_a, ...)`:
  * compute_adotoa (and, full-mapping, compute_background_EFT_functions) is called with the host `a` and the host
    `eft_cache` instead of `temp_a` / `temp_cache`, so temp_cache%adotoa stays 0 and the host cache is overwritten;
  * cP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2*eft_cache%EFTc/a (and LambdaP) uses the host cache and the
    host a, so it does not depend on temp_a: dfridr returns 0 and EFTcdotdot = -2*adotoa*EFTcdot;
  * full mapping only: the dummies gamma2PP, gamma3PP, gamma3PPP are never assigned (gfortran -Wuninitialized),
    so EFTGamma2PPP, EFTGamma3PPP, EFTGamma3PPPP are derivatives of uninitialised memory.
These feed compute_mass_factors (EFT_mu1, EFT_mu2), i.e. the EFT_mass_stability criterion.

Method: after camb.get_results the public EFTTestStability(a, k_max, EFTCAMB, par_cache, cache) is called through ctypes
(EFT_mass_stability is switched on after get_results so that compute_additional_derivs + compute_mass_factors run\nfor a > EFTCAMB_turn_on_time).
The returned EFTcdotdot / EFTGamma*PP / EFTLambdadotdot are compared with centred finite differences of the
first derivatives returned by the same routine at a*(1 +- eps):  GammaPP = d(GammaP)/da and, with the EFTCAMB
convention Xdot = a^2 dX/dtau (same as gpresdotm_t), Xdotdot = a*H*d(Xdot)/da - 2*H*Xdot.

Usage: /Users/marco/Software/miniforge3/bin/python additional-derivs-wrong-cache.py [TREE]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0"
)
PAR = REPO / "fortran/eftcamb_test/parameters"
CASES = [
    ("designer RPH", "3_RPH_KT_power_1.ini", {}),
    ("designer f(R)", "4_Des_F_R_wCDM_1.ini", {}),
    ("full mapping 5e quintessence", "5_quint_cubic.ini", {}),
]

CHILD = r"""
import sys, json, ctypes, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb
from camb.baseconfig import camblib
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.EFTCAMB.EFT_mass_stability = False   # run the model with its own stability settings minus the mass check
res = camb.get_results(pars)
P = res.Params
P.EFTCAMB.EFT_mass_stability = True     # then let EFTTestStability compute additional derivs + mass factors
f = camblib.__eftcamb_stability_MOD_eftteststability
f.restype = ctypes.c_int
def cache(a):
    c = eftcamb.EFTCAMB_timestep_cache()
    kmax = ctypes.c_double(10.0)
    ok = f(ctypes.byref(ctypes.c_double(a)), ctypes.byref(kmax), ctypes.byref(P.EFTCAMB),
           ctypes.byref(P.EFTCAMB_parameter_cache), ctypes.byref(c))
    return c, bool(ok)
out = []
for a in [0.1, 0.3, 0.6, 0.9]:
    c, ok = cache(a)
    eps = 1e-4
    cp, _ = cache(a * (1 + eps)); cm, _ = cache(a * (1 - eps))
    h = 2 * a * eps
    row = dict(a=a, stable=ok, adotoa=c.adotoa)
    # EFTCAMB convention (as for gpresdotm_t): Xdot = a^2 dX_phys/dtau, so Xdotdot = a*H*d(Xdot)/da - 2*H*Xdot
    row["EFTcdotdot"] = [c.EFTcdotdot, a * c.adotoa * (cp.EFTcdot - cm.EFTcdot) / h - 2 * c.adotoa * c.EFTcdot,
                         -2 * c.adotoa * c.EFTcdot]
    row["EFTLambdadotdot"] = [c.EFTLambdadotdot,
                              a * c.adotoa * (cp.EFTLambdadot - cm.EFTLambdadot) / h - 2 * c.adotoa * c.EFTLambdadot, 0.0]
    for g in ("1", "2", "3"):
        row[f"EFTGamma{g}PP"] = [getattr(c, f"EFTGamma{g}PP"), (getattr(cp, f"EFTGamma{g}P") - getattr(cm, f"EFTGamma{g}P")) / h, 0.0]
    row["EFTGamma2PPP"] = [c.EFTGamma2PPP, float("nan"), 0.0]
    row["EFT_mu1"] = [c.EFT_mu1, float("nan"), 0.0]
    row["EFT_mu2"] = [c.EFT_mu2, float("nan"), 0.0]
    out.append(row)
print("JSON" + json.dumps(out))
"""


def main():
    tmp = Path(tempfile.mkdtemp(prefix="addder_"))
    print(f"tree: {TREE}")
    for label, ini, extra in CASES:
        w = tmp / f"w_{ini}"
        lines = [f"DEFAULT({PAR / ini})", "feedback_level = 0", "EFTCAMB_write_background = F",
                 "get_tensor_cls = F", "do_lensing = F", f"output_root = {tmp / ini}"]
        lines += [f"{k} = {v}" for k, v in extra.items()]
        w.write_text("\n".join(lines) + "\n")
        r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w)], capture_output=True, text=True,
                           cwd=tmp, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
        js = [ln for ln in r.stdout.splitlines() if ln.startswith("JSON")]
        print(f"\n== {label}: {ini} (rc={r.returncode})")
        if not js:
            print(r.stdout[-1500:], r.stderr[-1500:])
            continue
        for row in json.loads(js[-1][4:]):
            print(f"  a = {row['a']:.2f}  stable={row['stable']}")
            for k in ("EFTcdotdot", "EFTLambdadotdot", "EFTGamma1PP", "EFTGamma2PP", "EFTGamma3PP"):
                code, fd, alt = row[k]
                rel = abs(code - fd) / max(abs(fd), 1e-300) if fd == fd else float("nan")
                extra_txt = f"   (-2*H*cdot = {alt:+.6e})" if k == "EFTcdotdot" else ""
                print(f"    {k:16s} code = {code:+.6e}   finite-diff = {fd:+.6e}   rel.err = {rel:.2e}{extra_txt}")
            print(f"    EFTGamma2PPP = {row['EFTGamma2PPP'][0]:+.6e}  EFT_mu1 = {row['EFT_mu1'][0]:+.6e}  "
                  f"EFT_mu2 = {row['EFT_mu2'][0]:+.6e}")
    print("\nexpected: code == finite-diff to ~1e-6 for every derivative")


if __name__ == "__main__":
    main()
