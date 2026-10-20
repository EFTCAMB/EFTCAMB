"""Adversarial verifier, finding t3-extgal-gr-point-background (independent check).

Independent of the repro (which compares harness npz files against an EFTflag=0 run by interpolation):
  * exact tracker root E^(2+S) = Omega_phi0 + Omega_tot(a) E^S solved here with brentq, Omega_tot(a) built from CAMB's
    own GR background densities (results.get_background_densities, incl. massive nu);
  * CAMB background H (results.hubble_parameter = dtauda) compared at exact z values with
      exact tracker H, and with the stated root cause prediction H0 sqrt(Omega_tot) (no DE) for a < 0.2;
  * comoving distance to z_star vs exact integral, for S = 0, 0.1, 0.5, 2 (massless nu, isolates the a < 0.2 branch);
  * get_eft_functions: NaN census at S = 0, and Hdotdot vs a finite difference of Hdot (S = 2, massless nu).
Each CAMB evaluation is a child process inside a harness CPU slot (OMP_NUM_THREADS=1).

Usage: python verify_extgal_gr_point.py [TREE] [OUTJSON]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
CHECKS = REPO / "fortran/eftcamb_test/audit/checks"
sys.path.insert(0, str(CHECKS))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work/verify/extgal_gr"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "get_scalar_cls": "F", "get_transfer": "F",
        "num_redshiftwindows": "0", "l_max_scalar": "200"}
C = 299792.458
Z_LIST = [0.5, 1.0, 3.0, 3.9, 3.99, 4.01, 4.1, 5.0, 9.0, 30.0, 100.0]
A_EFT = np.concatenate([[1e-10, 1e-6, 1e-3, 0.1, 0.19], np.linspace(0.3, 0.9, 601), [1.0]])

CHILD = r"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
cfg = json.loads(sys.argv[1])
sys.path.insert(0, cfg["checks"]); import harness
import numpy as np
with harness.slot(1):
    sys.path.insert(0, cfg["tree"]); import camb
    cp = camb.read_ini(cfg["ini"], no_validate=True)
    if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and cfg["eft"]:
        cp.EFTCAMB.initialize_parameters(cp, harness.parse_ini(cfg["ini"]), print_header=False)
    res = camb.get_background(cp)
    out = {"camb": camb.__file__}
    z = np.array(cfg["z"])
    out["H"] = list(map(float, res.hubble_parameter(z)))
    der = res.get_derived_params()
    out["derived"] = {k: float(v) for k, v in der.items()}
    out["chistar"] = float(res.comoving_radial_distance(der["zstar"]))
    out["H0"] = float(cp.H0)
    a = np.array(cfg["a_dens"])
    d = res.get_background_densities(a, vars=["cdm", "baryon", "photon", "neutrino", "nu", "de", "tot"])
    out["dens"] = {k: list(map(float, v)) for k, v in d.items()}
    if cfg["eft"]:
        fields, arr = cp.EFTCAMB.get_eft_functions(res, np.array(cfg["a_eft"]))
        out["eft"] = {f: list(map(float, arr[f])) for f in arr.dtype.names
                      if f in ("adotoa", "Hdot", "Hdotdot", "EFTc", "EFTLambda", "EFTcdot", "EFTLambdadot",
                               "EFTGamma1V", "EFTGamma1P", "EFTGamma2V", "EFTGamma2P")}
        out["eft_nan_fields"] = sorted(f for f in arr.dtype.names
                                       if arr[f].dtype.kind == "f" and not np.all(np.isfinite(arr[f])))
print("JSON" + json.dumps(out))
"""


def run_child(tree: str, name: str, ini: Path, ov: dict, eft: bool) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    wrapper = OUT / f"{name}.ini"
    lines = [f"DEFAULT({ini})", f"output_root = {OUT / name}", "EFTCAMB_write_background = F"]
    lines += [f"{k} = {v}" for k, v in {**FAST, **ov}.items()]
    wrapper.write_text("\n".join(lines) + "\n")
    a_dens = list(np.logspace(-4, 0, 4001))
    cfg = {"checks": str(CHECKS), "tree": str(harness.tree_path(tree)), "ini": str(wrapper), "eft": eft,
           "z": Z_LIST, "a_dens": a_dens, "a_eft": list(A_EFT)}
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, json.dumps(cfg)], capture_output=True, text=True, env=env,
                       cwd=OUT, timeout=900, preexec_fn=harness._nice)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if p.returncode != 0 or not line:
        return {"error": f"rc={p.returncode}", "tail": (p.stdout + p.stderr).splitlines()[-15:]}
    r = json.loads(line[-1][4:])
    r["a_dens"] = a_dens
    return r


def omega_tot_fun(gr: dict):
    """Omega_tot(a) (matter + radiation + massive nu, no DE) from CAMB GR densities (8 pi G rho a^4)."""
    a = np.array(gr["a_dens"])
    h0 = gr["H0"] / C
    s = sum(np.array(gr["dens"][k]) for k in ("cdm", "baryon", "photon", "neutrino", "nu"))
    om = s / (3 * h0**2 * a**4)
    la, lom = np.log(a), np.log(om)

    def f(x: float) -> float:
        return float(np.exp(np.interp(np.log(x), la, lom)))

    ophi0 = 1.0 - f(1.0)
    return f, ophi0, h0


def e_exact(a: float, s: float, f, ophi0: float) -> float:
    ot = f(a)
    g = lambda e: e ** (2 + s) - ophi0 - ot * e**s  # noqa: E731
    return brentq(g, np.sqrt(ot) * (1 - 1e-12), 10 * np.sqrt(ot + ophi0) + 10, xtol=1e-15, rtol=1e-15, maxiter=500)


def main() -> None:
    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    outjson = Path(sys.argv[2]) if len(sys.argv) > 2 else OUT / f"indep_{Path(str(tree)).name}.json"
    ext = {"EFTCAMB_skip_stability": "T"}
    jobs = {
        "gr_nu": (GR, {}, False),
        "gr_nu0": (GR, {"omnuh2": "0.0"}, False),
        "S0_nu": (INI, {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100", **ext}, True),
        "S0_nu0": (INI, {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100", "omnuh2": "0.0", **ext}, True),
        "S0.1_nu0": (INI, {"ExtendedGalileon_S": "0.1", "omnuh2": "0.0", **ext}, True),
        "S0.5_nu0": (INI, {"ExtendedGalileon_S": "0.5", "omnuh2": "0.0", **ext}, True),
        "S2_nu0": (INI, {"ExtendedGalileon_S": "2.0", "omnuh2": "0.0", **ext}, True),
    }
    with ThreadPoolExecutor(3) as ex:
        futs = {n: ex.submit(run_child, tree, n, ini, ov, eft) for n, (ini, ov, eft) in jobs.items()}
        res = {n: f.result() for n, f in futs.items()}
    report: dict = {"tree": str(harness.tree_path(tree))}
    for n, r in res.items():
        if "error" in r:
            print(n, "ERROR", r)
            report[n] = r
    funs = {"nu": omega_tot_fun(res["gr_nu"]), "nu0": omega_tot_fun(res["gr_nu0"])}
    print(f"tree: {report['tree']}")
    for n in ("S0_nu", "S0_nu0", "S0.1_nu0", "S0.5_nu0", "S2_nu0"):
        r = res[n]
        if "error" in r:
            continue
        s = float(n.split("_")[0][1:])
        f, ophi0, h0 = funs[n.split("_")[1]]
        rep = {"S": s, "Omega_phi0": ophi0}
        print(f"== {n}: S={s}, Omega_phi0={ophi0:.6f}")
        print("   z        H_code/H_exact-1   H_code/H_noDE-1 (root-cause prediction for a<0.2)")
        rows = []
        for z, hc in zip(Z_LIST, r["H"]):
            a = 1 / (1 + z)
            he = h0 * e_exact(a, s, f, ophi0) * C
            hn = h0 * np.sqrt(f(a)) * C
            rows.append((z, hc / he - 1, hc / hn - 1))
            print(f"   {z:7.2f}  {hc / he - 1:+.3e}         {hc / hn - 1:+.3e}")
        rep["H_rows"] = rows
        zs = r["derived"]["zstar"]
        ast = 1 / (1 + zs)
        chi_ex = quad(lambda la: np.exp(-la) / (h0 * e_exact(np.exp(la), s, f, ophi0)), np.log(ast), 0.0,
                      limit=400, epsabs=0, epsrel=1e-11, points=[np.log(0.2)])[0]
        rep["chistar_code"], rep["chistar_exact"] = r["chistar"], chi_ex
        print(f"   chi_star code/exact-1 = {r['chistar'] / chi_ex - 1:+.3e}  (code {r['chistar']:.3f} Mpc)")
        rep["thetastar"] = r["derived"]["thetastar"]
        if "gr" not in n:
            g = res["gr_nu" if n.endswith("_nu") else "gr_nu0"]
            if s == 0.0:
                rep["thetastar_vs_GR"] = r["derived"]["thetastar"] / g["derived"]["thetastar"] - 1
                rep["age_vs_GR"] = r["derived"]["age"] / g["derived"]["age"] - 1
                print(f"   S=0 vs EFTflag=0: thetastar ratio-1 = {rep['thetastar_vs_GR']:+.3e}, "
                      f"age ratio-1 = {rep['age_vs_GR']:+.3e}, chistar ratio-1 = {r['chistar'] / g['chistar'] - 1:+.3e}")
        e = r.get("eft", {})
        rep["eft_nan_fields"] = r.get("eft_nan_fields")
        print(f"   get_eft_functions non-finite fields: {r.get('eft_nan_fields')}")
        if e:
            aa = A_EFT
            ad = np.array(e["adotoa"])
            ex = np.array([a_ * h0 * e_exact(a_, s, f, ophi0) for a_ in aa])
            for av in (0.19, 0.3, 0.5, 1.0):
                i = int(np.argmin(abs(aa - av)))
                print(f"   eft adotoa/exact-1 at a={aa[i]:.3f}: {ad[i] / ex[i] - 1:+.3e}")
            # Hdotdot vs finite difference of Hdot on the uniform a grid 0.3..0.9 (conformal: d/dtau = a adotoa d/da)
            m = slice(5, 5 + 601)
            a_u, hd, hdd = aa[m], np.array(e["Hdot"])[m], np.array(e["Hdotdot"])[m]
            fd = a_u * ad[m] * np.gradient(hd, a_u, edge_order=2)
            # also Hdot vs finite difference of adotoa
            fd1 = a_u * ad[m] * np.gradient(ad[m], a_u, edge_order=2)
            k = slice(10, -10)
            if np.all(np.isfinite(hdd[k])):
                rel = np.max(abs(hdd[k] / fd[k] - 1))
                rel1 = np.max(abs(hd[k] / fd1[k] - 1))
                rep["Hdot_vs_fd_maxrel"], rep["Hdotdot_vs_fd_maxrel"] = float(rel1), float(rel)
                j = int(np.argmin(abs(a_u - 0.5)))
                print(f"   Hdot   vs a*adotoa*d(adotoa)/da: max rel diff (0.3<a<0.9) = {rel1:.2e}")
                print(f"   Hdotdot vs a*adotoa*d(Hdot)/da  : max rel diff (0.3<a<0.9) = {rel:.2e}"
                      f"  [a=0.5: code {hdd[j]:.4e}, FD {fd[j]:.4e}]")
            else:
                print("   Hdotdot non-finite:", hdd[k][:3])
            for fld in ("Hdotdot", "EFTcdot", "EFTLambdadot", "EFTGamma1P", "EFTc", "EFTLambda"):
                v = np.array(e[fld])
                print(f"   {fld:13s} a=0.5: {v[int(np.argmin(abs(aa - 0.5)))]:+.4e}   a=0.1: {v[3]:+.4e}")
        report[n] = rep
    outjson.write_text(json.dumps(report, indent=1, default=float))
    print("wrote", outjson)


if __name__ == "__main__":
    main()
