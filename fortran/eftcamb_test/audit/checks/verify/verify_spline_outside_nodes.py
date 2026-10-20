"""Adversarial verification of finding 'spline-derivatives-outside-nodes' (t4a-param).

Independent of the t4a repro (which compiles a Fortran driver and compares a single node layout):
 (a) attribution: the same constant w = -1 spline (EFTwDE = 9 / 10, null_value = -1) on node layouts that do or do
     not leave the range [x_1, x_N] during the run.  If the reported spectra difference is caused by the
     out-of-range branch, full coverage must give ~0 difference vs EFTwDE = 0, and the difference must come from the
     uncovered part of the a-range only.
 (b) physics consistency on the harness EFT-function grid: EFTLambdadot must equal d EFTLambda / d tau =
     a * adotoa * d EFTLambda / da (same for EFTcdot vs EFTc).  A spurious w' = null_value breaks this outside the nodes.

usage: python verify_spline_outside_nodes.py [run|analyse]   (default: both)
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "spline_outside"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"


def spline_keys(flag: int, nodes, null_value=-1.0, value=-1.0) -> dict:
    d = {"EFTwDE": str(flag), "EFTw_Spline_Pixels": str(len(nodes)), "EFTw_null_value": repr(null_value)}
    if flag == 10:
        d.update({"EFTw_Spline_Order": "5", "EFTw_Logx_Flag": "F", "EFTw_Logy_Flag": "F"})
    d.update({f"EFTwx{i + 1}": repr(float(x)) for i, x in enumerate(nodes)})
    d.update({f"EFTwv{i + 1}": repr(value) for i in range(len(nodes))})
    return d


FULL = np.logspace(-12, np.log10(2.0), 10)  # covers every a EFTCAMB can ask for (a_ini 1e-8 ... 1)
LATE_GAP = np.concatenate([np.logspace(-12, -2, 4), [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]])  # outside only a > 0.9
LATE_GAP2 = np.concatenate([np.logspace(-12, -2, 4), [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]])  # outside only a > 0.5
EARLY_GAP = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])  # flag 10: outside only a < 0.01 (x_final inclusive)
T4A = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9])

CASES = {
    "ref": {},
    "s9_t4a": spline_keys(9, T4A),
    "s10_t4a": spline_keys(10, T4A),
    "s9_full": spline_keys(9, FULL),
    "s10_full": spline_keys(10, FULL),
    "s10_lategap": spline_keys(10, LATE_GAP),
    "s10_lategap05": spline_keys(10, LATE_GAP2),
    "s10_earlygap": spline_keys(10, EARLY_GAP),
    "s9_endpoint1": spline_keys(9, np.logspace(-12, 0, 10)),  # last node exactly a=1 (cubic uses >=)
    # control: the bug is invisible when null_value = 0 (default) and the function is 0 outside
    # (w = 0 is not physical for DE, so instead use a constant w=-1 inside but check derivatives only via (b))
}


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    from concurrent.futures import ThreadPoolExecutor

    def one(item):
        name, ov = item
        return name, harness.run_model("A0", INI, OUT / f"{name}.npz", overrides=ov, threads=1)

    todo = {k: v for k, v in CASES.items() if not (OUT / f"{k}.npz").exists()}
    with ThreadPoolExecutor(4) as ex:
        for name, st in ex.map(one, todo.items()):
            print(f"{name:16s} {st['status']}  {st['elapsed']}s")


def dtau_numeric(a, f, adotoa):
    """a * adotoa * df/da with a central difference in ln a."""
    lna = np.log(a)
    dfdlna = np.gradient(f, lna)
    return adotoa * dfdlna  # d/dtau = adotoa * d/dln a


def consistency(d, field, dotfield, amin, amax):
    """EFTCAMB convention: X = a^2 x (conformal), Xdot = a^2 dx/dtau, so Xdot = a^2 d(X/a^2)/dtau."""
    a = d["eft_a_grid"]
    m = (a > amin) & (a < amax)
    num = a**2 * dtau_numeric(a, d[f"eft_{field}"] / a**2, d["eft_adotoa"])
    ana = d[f"eft_{dotfield}"]
    scale = np.max(np.abs(ana[m]))
    return float(np.max(np.abs(num[m] - ana[m])) / scale)


def analyse():
    ref = harness.load(OUT / "ref.npz")
    rows = {}
    print(f"{'case':16s} {'worstCl':>10s} {'TT':>10s} {'TE':>10s} {'P(k)':>10s} {'Pk@0.1':>10s} {'sigma8':>10s} "
          f"{'Ldot rel':>10s}")
    for name in CASES:
        if name == "ref":
            continue
        d = harness.load(OUT / f"{name}.npz")
        c = harness.compare(d, ref)
        k = d["pk_k"]
        j = np.argmin(abs(k - 0.1))
        pk01 = float(d["pk"][0, j] / ref["pk"][0, j] - 1)
        tt = c["cl"]["total"]["TT"]
        te = c["cl"]["total"]["TE"]
        rows[name] = {"worst_cl": c["worst_cl"], "TT": tt, "TE": te, "worst_pk": c["worst_pk"], "pk_k0.1": pk01,
                      "sigma8": c["sigma8"], "EFTLambdadot": c["eft"]["EFTLambdadot"],
                      "bit_identical": c["bit_identical"]}
        print(f"{name:16s} {c['worst_cl']:10.3e} {tt:10.3e} {te:10.3e} {c['worst_pk']:10.3e} {pk01:10.3e} "
              f"{c['sigma8']:10.3e} {c['eft']['EFTLambdadot']:10.3e}  bit_identical={c['bit_identical']}")
    print()
    print("(b) consistency: max |d/dtau X (numeric) - Xdot (code)| / max|Xdot|, on the window")
    windows = {"a in (0.02,0.85)": (0.02, 0.85), "a in (0.905,0.99)": (0.905, 0.99), "a=1 (edge)": (0.999, 1.0001), "a < 0.009": (1e-6, 0.009)}
    for name in ["ref", "s9_t4a", "s10_t4a", "s10_full", "s10_lategap", "s10_earlygap", "s9_endpoint1"]:
        d = ref if name == "ref" else harness.load(OUT / f"{name}.npz")
        parts = []
        for wl, (lo, hi) in windows.items():
            ld = consistency(d, "EFTLambda", "EFTLambdadot", lo, hi)
            cd = consistency(d, "EFTc", "EFTcdot", lo, hi)
            parts.append(f"{wl}: Ldot {ld:8.2e} cdot {cd:8.2e}")
            rows.setdefault(name, {})[f"consistency {wl}"] = {"EFTLambdadot": ld, "EFTcdot": cd}
        print(f"{name:14s} " + " | ".join(parts))
    # direct look at a = 1: the code's EFTLambdadot vs the w'-free expectation
    a = ref["eft_a_grid"]
    i1 = len(a) - 1
    s9 = harness.load(OUT / "s9_t4a.npz")
    print()
    print(f"a=1: EFTLambdadot ref={ref['eft_EFTLambdadot'][i1]:.6e}  s9_t4a={s9['eft_EFTLambdadot'][i1]:.6e}  "
          f"diff={s9['eft_EFTLambdadot'][i1] - ref['eft_EFTLambdadot'][i1]:.6e}  "
          f"expected diff if w'=-1: grhov*adotoa*a*(-1)={-ref['eft_grhov_t'][i1] * ref['eft_adotoa'][i1]:.6e}")
    print(f"a=1: EFTcdot      ref={ref['eft_EFTcdot'][i1]:.6e}  s9_t4a={s9['eft_EFTcdot'][i1]:.6e}  "
          f"diff={s9['eft_EFTcdot'][i1] - ref['eft_EFTcdot'][i1]:.6e}  "
          f"expected diff if w'=-1: 0.5*grhov*adotoa*a*(-1)={-0.5 * ref['eft_grhov_t'][i1] * ref['eft_adotoa'][i1]:.6e}")
    # decompose the a=1 EFTLambdadot change: explicit w' term + Hdotdot/Hdot change (DE pressure derivative uses w')
    O = ref["eft_EFTOmegaV"][i1]
    OP = ref["eft_EFTOmegaP"][i1]
    H = ref["eft_adotoa"][i1]
    dHdd = s9["eft_Hdotdot"][i1] - ref["eft_Hdotdot"][i1]
    dHd = s9["eft_Hdot"][i1] - ref["eft_Hdot"][i1]
    pred = -ref["eft_grhov_t"][i1] * H - 2 * O * (dHdd - H * dHd) - OP * (dHdd + 5 * H * dHd)
    print(f"a=1: dHdot={dHd:.3e} dHdotdot={dHdd:.3e} (Hdotdot ref {ref['eft_Hdotdot'][i1]:.3e}); predicted dLdot "
          f"from w'=-1 + dHdot/dHdotdot = {pred:.6e} vs observed "
          f"{s9['eft_EFTLambdadot'][i1] - ref['eft_EFTLambdadot'][i1]:.6e}")
    gq = s9["eft_gpresdotq"][i1] - ref["eft_gpresdotq"][i1]
    print(f"a=1: d gpresdotq = {gq:.6e}  (expected grhov*adotoa*a*w' = {-ref['eft_grhov_t'][i1] * H:.6e})")
    (OUT / "summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "both"
    if what in ("run", "both"):
        run()
    if what in ("analyse", "both"):
        analyse()
