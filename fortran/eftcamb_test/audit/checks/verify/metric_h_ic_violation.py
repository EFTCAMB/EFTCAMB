"""Independent check for finding metric-h-ic-violates-constraint (t5-robust).

Method (different from the repro, which compares spectra): use the instrumented tree work/verify/mh_instr (built by
metric_h_build_instr.py; inert without env vars, verified bitwise vs A0) whose derivs logs, for k ~ 1e-3..10 /Mpc,
k, tau, a, h'(used), h'_00 (00 constraint recomputed in derivs), h''(ii), etak, dgrho, dphi, dphi', ...
For EFTCAMB_evolve_metric_h = T:
  * at the first derivs call after the EFT turn-on (the IC point), measure v_IC = h'_used/h'_00 - 1, and compare with
    the prediction of the stated root cause, i.e. that the IC sets h' = 2kz_GR = (2 k etak + dgrho)/adotoa while the
    constraint is h'_00 = (4 k etak + 2 dgrho/Meff2)/(2 - alphaB_BS)/adotoa (delta_phi = 0 at the IC):
        v_pred = (2 k etak + dgrho) (2 - alphaB_BS) / (4 k etak + 2 dgrho/Meff2) - 1
    with Meff2, alphaB from the cached EFT functions at the turn-on (harness npz, eft_*).
  * the same with EFT_MH_ICFIX=1 (IC from the constraint, as in the proposed fix): v_IC must be ~1e-16.
  * the violation along a for both runs (max over a in [1e-6, 1]).

Usage: python metric_h_ic_violation.py [turn_on_time]
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

TREE = harness.AUDIT_ROOT / "work" / "verify" / "mh_instr"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "metric_h"
OUT.mkdir(parents=True, exist_ok=True)
INI = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
TON = sys.argv[1] if len(sys.argv) > 1 else "1e-8"
COLS = "k tau a hdot hc hdd etak dgrho dphi dphidot rsa tc nophot dgq etakdot sigma".split()
TARGETS = [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 1.0, 10.0]
CASES = {
    "T": {},
    "Ticfix": {"EFT_MH_ICFIX": "1"},
}


def run(item):
    name, env = item
    tag = f"icv_on{TON}_{name}"
    binf = OUT / f"{tag}.bin"
    npz = OUT / f"{tag}.npz"
    st = harness.run_model(TREE, INI, npz, overrides={"EFTCAMB_evolve_metric_h": "T", "EFTCAMB_turn_on_time": TON},
                           extra_env={**env, "EFT_MH_LOG": str(binf)}, threads=1, timeout=900)
    if st["status"] != "ok":
        raise RuntimeError(f"{name}: {st['status']} {st['log_tail'][-5:]}")
    raw = np.fromfile(binf, dtype=np.float64)
    raw = raw[: raw.size // 16 * 16].reshape(-1, 16)
    return name, ({c: raw[:, i] for i, c in enumerate(COLS)}, harness.load(npz), binf)


with ThreadPoolExecutor(2) as ex:
    res = dict(ex.map(run, CASES.items()))

summary = {"turn_on": float(TON)}
for name, (d, npz, binf) in res.items():
    ag = npz["eft_a_grid"]
    meff2 = float(np.interp(np.log(float(TON)), np.log(ag), npz["eft_Meff2"]))
    aB_bs = -2.0 * float(np.interp(np.log(float(TON)), np.log(ag), npz["eft_alphaB"]))
    print(f"case {name}: turn-on a={TON}, Meff2={meff2:.6f}, alphaB_BS={aB_bs:+.3e}")
    print(f"  {'k':>9s} {'a_IC':>9s} {'v_IC=h_used/h_00-1':>20s} {'v_pred(root cause)':>20s} "
          f"{'GR-term cancel':>15s} {'dec.max a>1e-6':>14s} {'dec.max a>1e-3':>14s}")
    rows = {}
    for kt in TARGETS:
        ks = np.unique(d["k"])
        kk = ks[np.argmin(abs(np.log(ks / kt)))]
        if abs(np.log(kk / kt)) > 0.03:
            continue
        m = d["k"] == kk
        o = np.argsort(d["tau"][m], kind="stable")
        x = {c: v[m][o] for c, v in d.items()}
        i0 = 0  # first derivs call after the switch (IC point)
        v_ic = x["hdot"][i0] / x["hc"][i0] - 1
        et, dg, k = x["etak"][i0], x["dgrho"][i0], x["k"][i0]
        v_pred = (2 * k * et + dg) * (2 - aB_bs) / (4 * k * et + 2 * dg / meff2) - 1
        cancel = abs(2 * k * et + dg) / max(abs(2 * k * et), abs(dg))  # how much h'_GR cancels (super-horizon)
        # violation per decade of a (h' oscillates and crosses zero
        # sub-horizon, so a pointwise ratio is meaningless there)
        dec = {}
        for lo in range(-7, 0):
            mm = (x["a"] >= 10.0**lo) & (x["a"] < 10.0 ** (lo + 1))
            if mm.any():
                # median over all derivs calls in the decade (robust to RK trial-stage states, which are logged
                # too and are off the solution by O(1) at early times), normalised by the decade's median |h'_00|
                dec[f"1e{lo}"] = float(np.median(np.abs(x["hdot"][mm] - x["hc"][mm])) / np.median(np.abs(x["hc"][mm])))
        late6 = max(v for kk_, v in dec.items() if int(kk_[2:]) >= -6)
        late3 = max(v for kk_, v in dec.items() if int(kk_[2:]) >= -3)
        rows[kt] = {"k": float(k), "a_ic": float(x["a"][i0]), "v_ic": float(v_ic), "v_pred": float(v_pred),
                    "cancel": float(cancel), "max_v_a_gt_1e-6": late6, "max_v_a_gt_1e-3": late3,
                    "dphi_ic": float(x["dphi"][i0]), "dphidot_ic": float(x["dphidot"][i0]), "decades": dec}
        print(f"  {k:9.4g} {x['a'][i0]:9.3e} {v_ic:+20.6e} {v_pred:+20.6e} {cancel:15.3e} {late6:14.2e} {late3:14.2e}")
        print("        per-decade median|h'-h'_00|/median|h'_00|: " + " ".join(f"{kk_}:{v:.1e}" for kk_, v in dec.items()))
    summary[name] = {"Meff2": meff2, "alphaB_BS": aB_bs, "rows": rows}
    binf.unlink()  # bulky log, analysed
(OUT / f"icv_on{TON}_summary.json").write_text(json.dumps(summary, indent=1))
