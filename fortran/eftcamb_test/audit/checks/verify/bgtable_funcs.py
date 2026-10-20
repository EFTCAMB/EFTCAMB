"""Verify t3-bg-table-resolution-gr-floor at the function level (independent of C_l / harness.compare).

For K-mouflage (eps2_0=-1e-8), ADE (f_ac=1e-8) and quintessence (exp. potential, lambda=3e-3) in A0:
  * evaluate EFTLambda(a) and H(z) on a fine grid near a ~ 0.5-1 (40 samples per table cell),
  * locate the table nodes x_i = x_ini + i h (h = (x_fin - x_ini)/(N-1)),
  * measure EFTLambda/a^2 (constant in the LCDM limit) at nodes and mid-cells -> expect relative bump cosh(h)-1,
  * measure dH/H vs EFTflag=0 (same cosmology) and compare with 0.5*Omega_DE(z)*linear-interp error of e^{2x},
  * repeat with 10x points.
Each configuration runs in its own subprocess under harness.slot(1), OMP_NUM_THREADS=1.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

P = harness.PARAMS_DIR
OUT = harness.AUDIT_ROOT / "work" / "verify" / "bgtable"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "2000", "num_redshiftwindows": "0",
        "get_transfer": "F"}

CONFIGS = {
    # name: (ini, overrides, x_ini, x_fin, N)
    "kmof_N1000": (P / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08"}, "log(1e-9)", "0", 1000),
    "kmof_N10000": (P / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08", "model_background_num_points": "10000"},
                    "log(1e-9)", "0", 10000),
    "gr_nu0": (P / "1_EFT_GR.ini", {"omnuh2": "0.0"}, None, None, None),
    "ade_N1500": (P / "5_ADE_1.ini", {"f_ac": "1e-08"}, "log(1e-10)", "0", 1500),
    "ade_N15000": (P / "5_ADE_1.ini", {"f_ac": "1e-08", "model_background_num_points": "15000"},
                   "log(1e-10)", "0", 15000),
    "quint_N1100": (P / "5_quint_exponential.ini", {"EFT_lambda": "0.003"}, "log(1e-8)", "log(1.1)", 1100),
    "quint_N11000": (P / "5_quint_exponential.ini", {"EFT_lambda": "0.003", "model_background_num_points": "11000"},
                     "log(1e-8)", "log(1.1)", 11000),
    "gr": (P / "1_EFT_GR.ini", {}, None, None, None),
}


def child(name: str) -> None:
    import numpy as np

    sys.path.insert(0, str(harness.TREES["A0"]))
    import camb

    ini, ov, xi, xf, n = CONFIGS[name]
    OUT.mkdir(parents=True, exist_ok=True)
    w = OUT / f"{name}.ini"
    w.write_text("\n".join([f"DEFAULT({ini})", f"output_root = {OUT / name}", "EFTCAMB_write_background = F"]
                           + [f"{k} = {v}" for k, v in {**FAST, **ov}.items()]) + "\n")
    cp = camb.read_ini(str(w), no_validate=True)
    res = camb.get_background(cp)
    out = {"name": name}
    # fine x grid for a in [0.45, 1]
    if n is not None:
        x_ini, x_fin = eval(xi, {"log": np.log}), eval(xf, {"log": np.log})
        h = (x_fin - x_ini) / (n - 1)
        i0 = int(np.ceil((np.log(0.45) - x_ini) / h))
        i1 = int(np.floor((0.0 - x_ini) / h))
        nodes = x_ini + h * np.arange(i0, i1 + 1)
        sub = np.linspace(0, 1, 41)[:-1]
        x = (nodes[:-1, None] + h * sub[None, :]).ravel()
        x = np.concatenate([x, [nodes[-1]]])
        a = np.exp(x)
        _, arr = cp.EFTCAMB.get_eft_functions(res, a)
        lam = np.asarray(arr["EFTLambda"])
        R = lam / a**2
        # values at nodes vs mid-cell
        Rn = R[::40]
        Rm = R[20::40]
        bump = Rm / (0.5 * (Rn[:-1] + Rn[1:])) - 1.0
        out.update(h=h, cosh_h_minus_1=float(np.cosh(h) - 1.0), lam_mid_bump_mean=float(np.mean(bump)),
                   lam_mid_bump_min=float(np.min(bump)), lam_mid_bump_max=float(np.max(bump)),
                   node_R_spread=float(np.ptp(Rn) / np.mean(np.abs(Rn))))
        # exact expected sawtooth for e^{2x} under linear interpolation
        mu = (x - x_ini) / h - np.floor((x - x_ini) / h + 1e-12)
        xl = x - mu * h
        lin_err = ((1 - mu) * np.exp(2 * xl) + mu * np.exp(2 * (xl + h))) / np.exp(2 * x) - 1.0
        np.savez(OUT / f"{name}_fine.npz", x=x, a=a, R=R, lin_err=lin_err, adotoa=np.asarray(arr["adotoa"]))
        z = 1 / a - 1
    else:
        z = None
    zf = np.linspace(0.0, 1.3, 1301)
    out["H_zf"] = list(map(float, res.hubble_parameter(zf)))
    out["Ode_zf"] = list(map(float, res.get_Omega("de", zf)))
    if z is not None:
        out["H_fine"] = list(map(float, res.hubble_parameter(z)))
        out["z_fine"] = list(map(float, z))
    (OUT / f"{name}.json").write_text(json.dumps(out))
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}))


def run_child(name: str) -> None:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        r = subprocess.run([sys.executable, __file__, "_child", name], env=env, capture_output=True, text=True,
                           timeout=600)
    if r.returncode != 0:
        print(name, "FAILED", r.returncode, r.stdout[-2000:], r.stderr[-3000:])


def main() -> None:
    import numpy as np
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(3) as ex:
        list(ex.map(run_child, CONFIGS))
    ld = {n: json.loads((OUT / f"{n}.json").read_text()) for n in CONFIGS if (OUT / f"{n}.json").exists()}
    zf = np.linspace(0.0, 1.3, 1301)
    print(f"{'config':14s} {'h':>8s} {'cosh(h)-1':>10s} {'Lam bump mid':>13s} {'max dH/H zf':>12s} "
          f"{'mean dH/H':>10s} {'max dH/H fine':>13s} {'pred max':>9s} {'corr(meas,pred)':>15s}")
    for n in CONFIGS:
        if n.startswith("gr") or n not in ld:
            continue
        ref = ld["gr_nu0" if n.startswith("kmof") else "gr"]
        d = ld[n]
        dH = np.array(d["H_zf"]) / np.array(ref["H_zf"]) - 1
        f = np.load(OUT / f"{n}_fine.npz")
        zfine = np.array(d["z_fine"])
        Href = np.interp(zfine, zf, np.array(ref["H_zf"]))
        # reference H on the fine z: interpolate the reference's smooth H (error ~(dz)^2 H'' ~ 1e-9)
        dHf = np.array(d["H_fine"]) / Href - 1
        Ode = np.interp(zfine, zf, np.array(ref["Ode_zf"]))
        pred = 0.5 * Ode * f["lin_err"]
        c = float(np.corrcoef(dHf, pred)[0, 1])
        print(f"{n:14s} {d['h']:8.5f} {d['cosh_h_minus_1']:10.3e} {d['lam_mid_bump_mean']:13.3e} "
              f"{np.max(np.abs(dH)):12.3e} {np.mean(dH):10.3e} {np.max(np.abs(dHf)):13.3e} {np.max(pred):9.3e} "
              f"{c:15.4f}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_child":
        child(sys.argv[2])
    else:
        main()
