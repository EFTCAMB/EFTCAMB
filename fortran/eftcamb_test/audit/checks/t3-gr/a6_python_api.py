"""t3-gr (a): EFTflag=0 Python-API paths changed by EFTCAMB (camb.py/model.py/results.py) and lensing sanity check.

python a6_python_api.py            -> runs worker on A0 and B0 in subprocesses and prints a comparison
python a6_python_api.py _worker TREE OUT.json
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def worker(tree: str, out: str) -> None:
    sys.path.insert(0, tree)
    import numpy as np

    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(tree).resolve())
    res: dict = {}
    base = dict(ombh2=0.0224, omch2=0.12, mnu=0.06, As=2.1e-9, ns=0.965, tau=0.055, lmax=2500)

    # 1) theta -> H0 solver (model.py: lower bracket raised to 1.1*100*sqrt(omega_m))
    for tag, kw in {
        "flat_thetastar": dict(thetastar=0.0104092, omk=0.0),
        "flat_cosmomc": dict(cosmomc_theta=0.0104109, omk=0.0),
        "closed_m01": dict(thetastar=0.0104092, omk=-0.1),
        "closed_m02": dict(thetastar=0.0104092, omk=-0.2),
        "flat_lowtheta": dict(thetastar=0.0098, omk=0.0),
    }.items():
        try:
            pars = camb.set_params(**base, **kw)
            r = camb.get_results(pars)
            cl = r.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]
            res[tag] = {"H0": pars.H0, "omega_de": r.get_Omega("de"), "TT1000": float(cl[1000, 0]),
                        "EE1000": float(cl[1000, 1]), "EFTflag": getattr(getattr(pars, "EFTCAMB", None), "EFTflag", None)}
        except Exception as e:  # noqa: BLE001
            res[tag] = {"error": repr(e)}

    # 2) matter power extrapolation beyond kmax (results.py: slope from 10-point average instead of last 2 points)
    pars = camb.set_params(H0=67.5, **base)
    pars.set_matter_power(redshifts=[0.0, 1.0], kmax=2.0)
    pars.NonLinear = camb.model.NonLinear_none
    r = camb.get_results(pars)
    kh = np.array([1.0, 1.9, 3.0, 10.0, 30.0, 100.0])
    for extrap in (True,):
        PK = r.get_matter_power_interpolator(nonlinear=False, hubble_units=True, k_hunit=True, extrap_kmax=200.0,
                                             log_interp=extrap)
        res["extrap_log"] = {"kh": kh.tolist(), "P_z0": PK.P(0.0, kh).tolist(), "P_z1": PK.P(1.0, kh).tolist()}
    try:
        PKn = r.get_matter_power_interpolator(nonlinear=False, hubble_units=True, k_hunit=True, extrap_kmax=200.0,
                                              log_interp=False)
        res["extrap_lin"] = {"P_z0": PKn.P(0.0, kh).tolist()}
    except Exception as e:  # noqa: BLE001
        res["extrap_lin"] = {"error": repr(e)}
    # the transfer k-grid spacing near kmax (non-uniform log spacing matters for the 10-point average slope)
    kk, zz, pk = r.get_linear_matter_power_spectrum(hubble_units=True, k_hunit=True)
    dl = np.diff(np.log(kk[-12:]))
    res["kgrid_last_dlogk"] = dl.tolist()

    # 3) lensing sanity check (lensing.f90: AmplitudeError commented out)
    for As in (1e-7, 1e-6):
        try:
            pars = camb.set_params(H0=67.5, **{**base, "As": As})
            r = camb.get_results(pars)
            cls = r.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)
            res[f"lens_As{As:g}"] = {"lensed_TT1000": float(cls["total"][1000, 0]),
                                     "unlensed_TT1000": float(cls["unlensed_total"][1000, 0]),
                                     "lensed_BB500": float(cls["total"][500, 2])}
        except Exception as e:  # noqa: BLE001
            res[f"lens_As{As:g}"] = {"error": repr(e)}

    # 4) set_params on an existing params object (camb.py: EFTCAMB re-initialised from the params dict each call)
    try:
        pars = camb.set_params(H0=67.5, **base)
        if hasattr(pars, "EFTCAMB"):
            pars.EFTCAMB.initialize_parameters(pars, {"EFTflag": 3, "DesignerEFTmodel": 1, "EFTwDE": 0,
                                                      "EFTB0": 0.01}, print_header=False)
            before = int(pars.EFTCAMB.EFTflag)
            camb.set_params(cp=pars, H0=68.0)
            res["set_params_existing"] = {"EFTflag_before": before, "EFTflag_after": int(pars.EFTCAMB.EFTflag)}
    except Exception as e:  # noqa: BLE001
        res["set_params_existing"] = {"error": repr(e)}

    Path(out).write_text(json.dumps(res, indent=1))


def main() -> None:
    from common import AUDIT, WORK

    outs = {}
    for t in ("A0", "B0"):
        out = WORK / f"a6_{t}.json"
        env = dict(os.environ, OMP_NUM_THREADS="2")
        env.pop("PYTHONPATH", None)
        p = subprocess.run([sys.executable, __file__, "_worker", str(AUDIT / "trees" / t), str(out)], env=env,
                           capture_output=True, text=True, cwd=WORK)
        if p.returncode:
            print(t, "failed", p.stdout[-2000:], p.stderr[-3000:])
        outs[t] = json.loads(out.read_text())
    a, b = outs["A0"], outs["B0"]
    for k in a:
        print(f"{k:22s} A0={json.dumps(a[k])[:300]}")
        print(f"{'':22s} B0={json.dumps(b.get(k))[:300]}")
    import numpy as np

    for key in ("extrap_log",):
        pa, pb = np.array(a[key]["P_z0"]), np.array(b[key]["P_z0"])
        print(key, "kh", a[key]["kh"], "A0/B0-1 =", (pa / pb - 1).round(6).tolist())


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(sys.argv[2], sys.argv[3])
    else:
        main()
