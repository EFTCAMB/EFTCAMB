"""Track 0: value round trips through the ctypes structures, for designer f(R), Horndeski hdsk_TG and K-mouflage.

For each model (run in its own subprocess, camb imported from the chosen tree):
  * parameter cache read through ctypes vs CAMBdata fields read through ctypes and vs analytic values,
  * get_eft_functions: echo of a and tau, adotoa vs a H(z)/c from results.hubble_parameter, species densities vs
    the analytic a-scaling of the parameter-cache densities, grhom_t and grhoa2 sums, and finite-difference checks of
    the derivative fields (catches swapped/misplaced fields),
  * designer f(R): B0 reconstructed from EFTOmegaV/EFTOmegaP/adotoa/Hdot at a=1 vs the input EFTB0,
  * get_scale_evolution / get_evolution: array shape and [k][time] order, echo of k, a and tau, consistency of the
    background fields with get_eft_functions, behaviour for decreasing a / unsorted eta (np.sort inside get_evolution).
Results are written as JSON to the work directory.

Usage: python roundtrip.py [--tree PATH] [--jobs 3]           (driver)
       python roundtrip.py --worker MODEL_INI OUT_JSON --tree PATH
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
WORK = ROOT / "work" / "t0-layout" / "roundtrip"
MODELS = ["4_Des_F_R_LCDM_1", "5_hdsk_TG_1", "5_Kmouflage_1"]


def rel(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    den = np.maximum(np.abs(b), 1e-300)
    return float(np.max(np.abs(a - b) / den))


def worker(model_ini: Path, out_json: Path, tree: Path) -> None:
    sys.path.insert(0, str(tree))
    import camb
    from camb import constants

    assert Path(camb.__file__).resolve().is_relative_to(tree.resolve())
    out_json.parent.mkdir(parents=True, exist_ok=True)
    wrapper = out_json.with_suffix(".ini")
    wrapper.write_text(
        f"DEFAULT({model_ini})\noutput_root = {out_json.with_suffix('')}\nEFTCAMB_write_background = F\n"
        "feedback_level = 0\n"
    )
    checks: dict[str, dict] = {}

    def check(name, value, thr, extra=None):
        checks[name] = {"value": float(value), "threshold": thr, "pass": bool(value <= thr)}
        if extra:
            checks[name].update(extra)

    cp = camb.read_ini(str(wrapper), no_validate=True)
    res = camb.get_background(cp)
    c_km = constants.c / 1000.0
    # reference values from the CAMBdata copy of the parameters (CAMB may change H0, num_nu_massive in SetParams)
    rp = res.Params
    H0, h = rp.H0, rp.H0 / 100.0
    checks["_input_vs_used"] = {"input_H0": cp.H0, "used_H0": rp.H0, "input_num_nu_massive": cp.num_nu_massive,
                                "used_num_nu_massive": rp.num_nu_massive}

    # --- 1) parameter cache through ctypes vs CAMBdata through ctypes and analytic values
    pc = res.Params.EFTCAMB_parameter_cache
    # the parameter cache holds the model's H0 (Horndeski overwrites it with its own background value)
    _, A1 = cp.EFTCAMB.get_eft_functions(res, np.array([1.0]))
    checks["_h0"] = {"pc.h0": pc.h0, "pc.h0_Mpc": pc.h0_Mpc, "CP.H0": H0, "adotoa(a=1)": float(A1["adotoa"][0]),
                     "hubble_parameter(0)": float(res.hubble_parameter(0.0))}
    check("pc.h0_Mpc == adotoa(a=1)", rel(pc.h0_Mpc, A1["adotoa"][0]), 1e-12)
    check("pc.h0 == hubble_parameter(z=0) [km/s/Mpc, not H0/100 as documented]", rel(pc.h0, res.hubble_parameter(0.0)), 1e-6)
    for fld in ("grhog", "grhornomass", "grhoc", "grhob", "grhov", "grhok"):
        a_, b_ = getattr(pc, fld), getattr(res, fld)
        check(f"pc.{fld} == CAMBdata.{fld}", abs(a_ - b_) / max(abs(b_), 1e-300) if b_ else abs(a_), 1e-14)
    grhom = 3 * (1e5 / constants.c) ** 2
    check("pc.grhoc == 3 (1e5/c)^2 omch2", rel(pc.grhoc, grhom * rp.omch2), 1e-12)
    check("pc.grhob == 3 (1e5/c)^2 ombh2", rel(pc.grhob, grhom * rp.ombh2), 1e-12)
    check("pc.omegac == omch2/h^2", rel(pc.omegac, rp.omch2 / h**2), 1e-14)
    check("pc.omegab == ombh2/h^2", rel(pc.omegab, rp.ombh2 / h**2), 1e-14)
    check("pc.omegan == omnuh2/h^2", abs(pc.omegan - rp.omnuh2 / h**2), 1e-14)
    check("pc.omegav == CAMBdata.omega_de", rel(pc.omegav, res.omega_de), 1e-14)
    check("pc.omegak == omk", abs(pc.omegak - rp.omk), 0.0)
    check("pc.omegag == grhog/grhocrit", rel(pc.omegag, res.grhog / res.grhocrit), 1e-14)
    check("pc.Num_Nu_Massive == num_nu_massive", abs(pc.Num_Nu_Massive - rp.num_nu_massive), 0)
    check("pc.Nu_mass_eigenstates == nu_mass_eigenstates", abs(pc.Nu_mass_eigenstates - rp.nu_mass_eigenstates), 0)
    gm, nm = np.array(pc.grhormass[:]), np.array(pc.nu_masses[:])
    check("pc.grhormass (AllocatableArrayDouble) == CAMBdata.grhormass", float(np.max(np.abs(gm - np.array(res.grhormass[:])))), 0.0,
          {"size": len(gm), "values": gm.tolist()})
    check("pc.nu_masses (AllocatableArrayDouble) == CAMBdata.nu_masses", float(np.max(np.abs(nm - np.array(res.nu_masses[:])))), 0.0,
          {"values": nm.tolist()})

    # --- 2) get_eft_functions
    a = np.logspace(-4, 0, 121)
    z = 1 / a - 1
    fields, arr = cp.EFTCAMB.get_eft_functions(res, a)
    names = list(arr.dtype.names)
    check("eft fields list == timestep cache fields", 0 if names == fields else 1, 0, {"n_fields": len(names)})
    check("eft itemsize == 1272", abs(arr.dtype.itemsize - 1272), 0)
    check("eft a echo", float(np.max(np.abs(arr["a"] - a))), 0.0)
    check("eft tau echo == conformal_time(z)", rel(arr["tau"], res.conformal_time(z)), 0.0)
    adot_ref = a * res.hubble_parameter(z) / c_km
    check("eft adotoa vs a H(z)/c", rel(arr["adotoa"], adot_ref), 1e-6)
    check("eft grhoc_t == grhoc/a", rel(arr["grhoc_t"], pc.grhoc / a), 1e-14)
    check("eft grhob_t == grhob/a", rel(arr["grhob_t"], pc.grhob / a), 1e-14)
    check("eft grhog_t == grhog/a^2", rel(arr["grhog_t"], pc.grhog / a**2), 1e-14)
    check("eft grhor_t == grhornomass/a^2", rel(arr["grhor_t"], pc.grhornomass / a**2), 1e-14)
    gsum = arr["grhob_t"] + arr["grhoc_t"] + arr["grhog_t"] + arr["grhor_t"] + arr["grhonu_tot"]
    check("eft grhom_t == b+c+g+r+nu", rel(arr["grhom_t"], gsum), 1e-13)
    check("eft grhoa2 == grhom_t/a^4", rel(arr["grhoa2"], arr["grhom_t"] / a**4), 1e-13)
    check("eft k == 0 (not set by get_eft_functions)", float(np.max(np.abs(arr["k"]))), 0.0)
    # neutrino density vs CAMB massive neutrino background: grhonu_tot = sum grhormass/a^2 * rho(a m)
    if pc.Num_Nu_Massive:
        from camb import model as _m  # noqa: F401
        check("eft grhonu_tot > 0 with massive nu", 0.0 if np.all(arr["grhonu_tot"] > 0) else 1.0, 0.0)

    # finite-difference derivative checks (field ordering)
    a0s = np.array([0.05, 0.2, 0.5, 0.9])
    hstep = 1e-4
    ap, am = a0s * (1 + hstep), a0s * (1 - hstep)
    _, A0 = cp.EFTCAMB.get_eft_functions(res, a0s)
    _, Ap = cp.EFTCAMB.get_eft_functions(res, ap)
    _, Am = cp.EFTCAMB.get_eft_functions(res, am)
    da = ap - am
    dtau_da = 1 / (a0s * A0["adotoa"])  # d tau / d a = 1/(a adotoa)
    pairs_a = [("EFTOmegaV", "EFTOmegaP"), ("EFTOmegaP", "EFTOmegaPP"), ("EFTOmegaPP", "EFTOmegaPPP"),
               ("EFTGamma1V", "EFTGamma1P"), ("EFTGamma2V", "EFTGamma2P"), ("EFTGamma3V", "EFTGamma3P"),
               ("EFTGamma4V", "EFTGamma4P"), ("EFTGamma5V", "EFTGamma5P"), ("EFTGamma6V", "EFTGamma6P")]
    pairs_tau = [("adotoa", "Hdot"), ("Hdot", "Hdotdot"), ("EFTc", "EFTcdot"),
                 ("EFTLambda", "EFTLambdadot"), ("alphaB", "alphaBdot"), ("alphaT", "alphaTdot"),
                 ("alphaK", "alphaKdot"), ("alphaM", "alphaMdot"), ("phi_scf", "dphi"), ("grhoq", "grhodotq"),
                 ("gpresq", "gpresdotq"), ("gpresm_t", "gpresdotm_t"), ("gpresdotm_t", "gpresdotdotm_t"),
                 ("grhonu_tot", "grhonudot_tot"), ("gpinu_tot", "gpinudot_tot"), ("gpinudot_tot", "gpinudotdot_tot")]
    a2conv = {"EFTc", "EFTLambda", "grhoq", "gpresq", "gpresm_t", "gpresdotm_t", "grhonu_tot", "gpinu_tot",
              "gpinudot_tot"}
    fd = {}
    for f, df in pairs_a + pairs_tau:
        num = (Ap[f] - Am[f]) / da
        if (f, df) in pairs_tau:
            num = num / dtau_da
        ana = A0[df]
        if f in a2conv:
            ana = ana + 2 * A0["adotoa"] * A0[f]
        scale = max(np.max(np.abs(ana)), np.max(np.abs(num)))
        if scale == 0:
            fd[f"d{f}"] = {"skipped": "identically zero"}
            continue
        if (f, df) in pairs_a and np.max(np.abs(A0[f])) < 1e-10:
            # dimensionless EFT function at numerical-noise level (e.g. Horndeski TG gamma1, gamma2 ~ 1e-13)
            fd[f"d{f}"] = {"skipped": "noise level", "max_abs_value": float(np.max(np.abs(A0[f]))),
                           "num": num.tolist(), "ana": ana.tolist()}
            continue
        err = float(np.max(np.abs(num - ana)) / scale)
        fd[f"d{f}"] = {"field": df, "max_rel_err": err, "num": num.tolist(), "ana": ana.tolist()}
    worst = max([v["max_rel_err"] for v in fd.values() if "max_rel_err" in v] or [0.0])
    worst_name = max((k for k in fd if "max_rel_err" in fd[k]), key=lambda k: fd[k]["max_rel_err"], default="")
    check("finite-difference derivative fields (worst)", worst, 2e-2, {"worst": worst_name})
    checks["_fd_detail"] = fd
    # Hdotdotdot separately (finding fullmap-hdotdotdot-wrong): numerical d(Hdotdot)/dtau vs the field
    num3 = (Ap["Hdotdot"] - Am["Hdotdot"]) / da / dtau_da
    check("Hdotdotdot vs d(Hdotdot)/dtau", float(np.max(np.abs(num3 - A0["Hdotdotdot"])) / np.max(np.abs(num3))), 2e-2,
          {"num": num3.tolist(), "field": A0["Hdotdotdot"].tolist()})

    # --- 3) designer f(R): B0 from the EFT functions
    ini_vals = {}
    for line in model_ini.read_text().splitlines():
        if "=" in line and not line.strip().startswith("#"):
            k_, v_ = line.split("=", 1)
            ini_vals[k_.strip()] = v_.strip()
    if ini_vals.get("EFTflag") == "3" and ini_vals.get("DesignerEFTmodel") == "1":
        _, A1 = cp.EFTCAMB.get_eft_functions(res, np.array([1.0]))
        om, omp, H, Hd = A1["EFTOmegaV"][0], A1["EFTOmegaP"][0], A1["adotoa"][0], A1["Hdot"][0]
        B0 = omp / (1 + om) * H**2 / (Hd - H**2)  # B = (dOmega/dlna)/(1+Omega) * H/(dH/dlna), a = 1
        B0_in = float(ini_vals["EFTB0"])
        check("designer f(R) B0 from EFT functions vs EFTB0", abs(B0 - B0_in) / abs(B0_in), 1e-2,
              {"B0_reconstructed": float(B0), "EFTB0": B0_in, "OmegaV(1)": float(om), "OmegaP(1)": float(omp)})

    # --- 4) get_scale_evolution / get_evolution
    ks = [0.2, 0.05]  # deliberately non-increasing
    a_ev = np.logspace(-3, 0, 25)
    fl, ev = cp.EFTCAMB.get_scale_evolution(res, ks, a_ev)
    check("evolution shape == (nk, na)", 0 if ev.shape == (len(ks), len(a_ev)) else 1, 0, {"shape": list(ev.shape)})
    check("evolution k echo [k][time]", float(np.max(np.abs(ev["k"] - np.array(ks)[:, None]))), 0.0)
    eta_ev = res.conformal_time(1 / a_ev - 1)
    check("evolution tau echo (increasing a)", rel(ev["tau"], np.broadcast_to(eta_ev, ev.shape)), 1e-12)
    check("evolution a echo (increasing a)", rel(ev["a"], np.broadcast_to(a_ev, ev.shape)), 1e-4)
    _, ef = cp.EFTCAMB.get_eft_functions(res, a_ev)
    m = a_ev >= max(cp.EFTCAMB.EFTCAMB_pert_turn_on, 1e-3)
    check("evolution adotoa == get_eft_functions adotoa", rel(ev["adotoa"][:, m], np.broadcast_to(ef["adotoa"][m], ev[:, m].shape)), 1e-4)
    check("evolution EFTOmegaV == get_eft_functions EFTOmegaV",
          float(np.max(np.abs(ev["EFTOmegaV"][:, m] - ef["EFTOmegaV"][m]))) / max(np.max(np.abs(ef["EFTOmegaV"][m])), 1e-300),
          1e-3)
    # decreasing a: the caller's order is lost because get_evolution sorts eta
    fl, evr = cp.EFTCAMB.get_scale_evolution(res, ks, a_ev[::-1])
    check("evolution a echo with decreasing a (caller order)", rel(evr["a"][0], a_ev[::-1]), 1e-4,
          {"returned_a_first_last": [float(evr["a"][0, 0]), float(evr["a"][0, -1])],
           "requested_a_first_last": [float(a_ev[::-1][0]), float(a_ev[::-1][-1])]})
    rng = np.random.default_rng(1)
    perm = rng.permutation(len(eta_ev))
    fl, evp = cp.EFTCAMB.get_evolution(res, ks, eta_ev[perm])
    check("get_evolution tau echo with unsorted eta (caller order)", rel(evp["tau"][0], eta_ev[perm]), 1e-12)
    # scalar inputs
    fl, evs = cp.EFTCAMB.get_evolution(res, 0.1, float(eta_ev[-1]))
    check("get_evolution scalar inputs shape (1,1)", 0 if evs.shape == (1, 1) else 1, 0)

    # --- 5) TEFTCAMB python <-> fortran: python fields after read_ini vs ini values
    flag_keys = {
        "EFTflag": "EFTflag", "DesignerEFTmodel": "DesignerEFTmodel", "FullMappingEFTmodel": "FullMappingEFTmodel",
        "EFT_ghost_math_stability": "EFT_ghost_math_stability", "EFT_mass_math_stability": "EFT_mass_math_stability",
        "EFT_ghost_stability": "EFT_ghost_stability", "EFT_gradient_stability": "EFT_gradient_stability",
        "EFT_mass_stability": "EFT_mass_stability", "EFT_mass_stability_rate": "EFT_mass_stability_rate",
        "EFT_additional_priors": "EFT_additional_priors", "EFTCAMB_stability_time": "EFTCAMB_stability_time",
        "EFTCAMB_stability_threshold": "EFTCAMB_stability_threshold",
        "EFTCAMB_evolve_delta_phi": "EFTCAMB_evolve_delta_phi", "EFTCAMB_evolve_metric_h": "EFTCAMB_evolve_metric_h",
        "EFTCAMB_skip_stability": "EFTCAMB_skip_stability", "EFTCAMB_skip_RGR": "EFTCAMB_skip_RGR",
    }
    bad = {}
    for key, attr in flag_keys.items():
        if key in ini_vals:
            v = ini_vals[key]
            pyv = getattr(cp.EFTCAMB, attr)
            exp = {"T": True, "F": False}.get(v)
            if exp is None:
                exp = float(v)
            if (isinstance(exp, bool) and pyv is not exp) or (not isinstance(exp, bool) and float(pyv) != exp):
                bad[key] = [str(pyv), v]
    check("TEFTCAMB python fields == ini values after read_ini", len(bad), 0, {"bad": bad})

    out_json.write_text(json.dumps(checks, indent=1))
    print(json.dumps({k: (v["value"], v["pass"]) for k, v in checks.items() if not k.startswith("_")}, indent=1))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", default=str(ROOT / "trees" / "A0"))
    p.add_argument("--worker", nargs=2)
    p.add_argument("--models", default=",".join(MODELS))
    p.add_argument("--tag", default="A0")
    a = p.parse_args()
    if a.worker:
        worker(Path(a.worker[0]).resolve(), Path(a.worker[1]), Path(a.tree))
        return
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    WORK.mkdir(parents=True, exist_ok=True)
    procs = []
    for mname in a.models.split(","):
        out = WORK / f"{a.tag}_{mname}.json"
        log = open(WORK / f"{a.tag}_{mname}.log", "w")
        procs.append((mname, out, subprocess.Popen(
            [sys.executable, __file__, "--tree", a.tree, "--worker", str(PARAMS / f"{mname}.ini"), str(out)],
            cwd=WORK, env=env, stdout=log, stderr=subprocess.STDOUT)))
    for mname, out, pr in procs:
        rc = pr.wait()
        print(f"== {mname}: rc={rc}")
        if rc == 0:
            for k, v in json.loads(out.read_text()).items():
                if not k.startswith("_"):
                    print(f"   {'PASS' if v['pass'] else 'FAIL'} {k}: {v['value']:.3e} (thr {v['threshold']})")


if __name__ == "__main__":
    main()
