"""Independent check of finding stability-gpinu-tot-tripled (verifier).

Does NOT use the python bridge (get_eft_functions) as reference. Instead, using only the stability cache filled by
__eftcamb_stability_MOD_eftteststability (i.e. EFTStabilityComputation):

  1. internal consistency: gpresm_t is assembled after the FIRST neutrino loop, so
     (gpresm_t - (grhor_t + grhog_t)/3) is the 1x massive-nu pressure; compare with the final gpinu_tot in the cache.
  2. physics: at a = 1e-6 massive neutrinos are ultra-relativistic, so w_nu = gpinu_tot / grhonu_tot must be ~1/3.
  3. finite differences in ln a of the stability-cache Hdot (computed from gpresm_t, so unaffected) give
     Hdotdot_FD = adotoa dHdot/dlna and Hdotdotdot_FD = adotoa d(adotoa dHdot/dlna)/dlna; compare with the cache values.

Usage: python gpinu_fd_indep.py OUT.json TREE [TREE ...]
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
AUD = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/inis")
MODELS = [
    PAR / "2_PEFT_gamma3_const_1.ini",
    PAR / "2_PEFT_wDE_wCDM_1.ini",
    PAR / "4_Des_F_R_LCDM_1.ini",
    PAR / "4_Des_mc5e_wCDM_1.ini",
    PAR / "5_Horava_eta_1.ini",
    PAR / "5_hdsk_TG_1.ini",  # control: full-mapping generic H derivs use gpresm_t only
    PAR / "5_quint_quadratic.ini",  # control
    AUD / "007p8-3p1.ini",
    AUD / "007p7-shiftsym.ini",
]

WORKER = r"""
import ctypes, json, sys
from ctypes import byref, c_double, c_int
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
from camb.baseconfig import camblib as lib
from camb.eftcamb import EFTCAMB_timestep_cache
cp = camb.read_ini(sys.argv[2], no_validate=True)
cp.EFTCAMB.EFTCAMB_skip_stability = True
res = camb.get_background(cp)
P = res.Params
test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
buf = (ctypes.c_char * (1 << 18))()
names = ["a", "adotoa", "Hdot", "Hdotdot", "Hdotdotdot", "gpinu_tot", "grhonu_tot", "gpresm_t", "grhor_t", "grhog_t"]
off = {n: getattr(EFTCAMB_timestep_cache, n).offset for n in names}
def cache(x):
    test(byref(c_double(x)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
    raw = bytes(buf[: max(off.values()) + 8])
    return {n: float(np.frombuffer(raw, dtype=np.float64, count=1, offset=off[n])[0]) for n in names}
h = 2e-3
out = {"camb": camb.__file__, "num_nu_massive": int(P.num_nu_massive), "omnuh2": float(P.omnuh2), "rows": []}
for a0 in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5]:
    cs = [cache(a0 * np.exp(j * h)) for j in (-2, -1, 0, 1, 2)]
    Hd = np.array([c["Hdot"] for c in cs]); ad = np.array([c["adotoa"] for c in cs])
    d1 = lambda f: (-f[4] + 8 * f[3] - 8 * f[1] + f[0]) / (12 * h)
    d2 = lambda f: (-f[4] + 16 * f[3] - 30 * f[2] + 16 * f[1] - f[0]) / (12 * h * h)
    c0 = cs[2]
    hdd_fd = c0["adotoa"] * d1(Hd)
    hddd_fd = c0["adotoa"] * (d1(ad) * d1(Hd) + c0["adotoa"] * d2(Hd))
    p1x = c0["gpresm_t"] - (c0["grhor_t"] + c0["grhog_t"]) / 3.0
    out["rows"].append(dict(a=a0, cache=c0, Hdotdot_fd=hdd_fd, Hdotdotdot_fd=hddd_fd, gpinu_1x_from_gpresm=p1x))
print("RESULT " + json.dumps(out))
"""


def run(tree, ini):
    with tempfile.TemporaryDirectory() as d:
        w = Path(d) / "m.ini"
        w.write_text(f"DEFAULT({ini})\noutput_root = {d}/out\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        with harness.slot(1):
            p = subprocess.run([sys.executable, "-c", WORKER, tree, str(w)], capture_output=True, text=True, cwd=d,
                               env=env, preexec_fn=lambda: os.nice(10), timeout=900)
        res = [json.loads(line[7:]) for line in p.stdout.splitlines() if line.startswith("RESULT ")]
        if not res:
            return {"error": (p.stdout[-1500:] + p.stderr[-1500:])}
        return res[0]


def main():
    out_json = sys.argv[1]
    trees = sys.argv[2:]
    allres = {}
    for tree in trees:
        for ini in MODELS:
            if not ini.exists():
                continue
            r = run(tree, str(ini))
            allres[f"{Path(tree).name}:{ini.stem}"] = r
            print(f"== tree {Path(tree).name}  model {ini.stem}")
            if "error" in r:
                print("  ERROR", r["error"][-600:])
                continue
            print("      a     gpinu_cache/gpinu_1x  w_nu=P/rho   (Hdd_cache-Hdd_FD)/|FD|  (Hddd_cache-Hddd_FD)/|FD|")
            for row in r["rows"]:
                c = row["cache"]
                ratio = c["gpinu_tot"] / row["gpinu_1x_from_gpresm"] if row["gpinu_1x_from_gpresm"] else float("nan")
                w = c["gpinu_tot"] / c["grhonu_tot"] if c["grhonu_tot"] else float("nan")
                e2 = (c["Hdotdot"] - row["Hdotdot_fd"]) / abs(row["Hdotdot_fd"])
                e3 = (c["Hdotdotdot"] - row["Hdotdotdot_fd"]) / abs(row["Hdotdotdot_fd"]) if c["Hdotdotdot"] else float(
                    "nan")
                print(f"   {row['a']:8.1e}   {ratio:9.5f}          {w:8.5f}      {e2: 11.3e}             {e3: 11.3e}")
    Path(out_json).write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()
