"""Verifier (t3-eft-turnon-before-tauminn-lmax): evaluate the refusal and two candidate behaviours built by
tauminn_build_fix.py (mode 1 = earlier thermal-history start [proposed fix], mode 2 = guard removed).

Checks: (a) modes 1/2 are bit-identical to A0 where A0 already runs; (b) where A0 refuses, modes 1/2 run and agree
with a nearby passing configuration to the same level as the GR (base_params) baseline for the same settings change.

Usage: python tauminn_fix_eval.py [run|report]
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "tauminn" / "fixeval"
FIX = harness.AUDIT_ROOT / "work" / "verify" / "tauminn_fix"
P = harness.PARAMS_DIR
NOTENS = {"get_tensor_cls": "F"}
NOWIN = {"num_redshiftwindows": "0"}
TONLY = {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_vector_cls": "F", "do_lensing": "F",
         "num_redshiftwindows": "0"}

RUNS = {
    # name: (tree, model, overrides, mode)
    "fix0_ADE_default": (FIX, "5_ADE_1", {}, "0"),
    "fix1_ADE_default": (FIX, "5_ADE_1", {}, "1"),
    "fix2_ADE_default": (FIX, "5_ADE_1", {}, "2"),
    "fix0_ADE_l2000": (FIX, "5_ADE_1", {**NOTENS, "l_max_scalar": "2000"}, "0"),
    "fix1_ADE_l2000": (FIX, "5_ADE_1", {**NOTENS, "l_max_scalar": "2000"}, "1"),
    "fix2_ADE_l2000": (FIX, "5_ADE_1", {**NOTENS, "l_max_scalar": "2000"}, "2"),
    "A0_ADE_l2500": ("A0", "5_ADE_1", {**NOTENS, "l_max_scalar": "2500"}, None),
    "A0_GR_l2000": ("A0", "base_params", {**NOTENS, "l_max_scalar": "2000"}, None),
    "A0_GR_l2500": ("A0", "base_params", {**NOTENS, "l_max_scalar": "2500"}, None),
    "fix1_ADE_nowin": (FIX, "5_ADE_1", {**NOTENS, **NOWIN}, "1"),
    "fix2_ADE_nowin": (FIX, "5_ADE_1", {**NOTENS, **NOWIN}, "2"),
    "A0_ADE_win_notens": ("A0", "5_ADE_1", NOTENS, None),
    "A0_GR_nowin": ("A0", "base_params", {**NOTENS, **NOWIN}, None),
    "A0_GR_win_notens": ("A0", "base_params", NOTENS, None),
    "fix0_quint_default": (FIX, "5_quint_cubic", {}, "0"),
    "fix1_quint_default": (FIX, "5_quint_cubic", {}, "1"),
    "fix1_quint_tonly_k0.1": (FIX, "5_quint_cubic", {**TONLY, "transfer_kmax": "0.1"}, "1"),
    "fix2_quint_tonly_k0.1": (FIX, "5_quint_cubic", {**TONLY, "transfer_kmax": "0.1"}, "2"),
    "A0_quint_tonly_k2": ("A0", "5_quint_cubic", TONLY, None),
    "A0_GR_tonly_k0.1": ("A0", "base_params", {**TONLY, "transfer_kmax": "0.1"}, None),
    "A0_GR_tonly_k2": ("A0", "base_params", TONLY, None),
}


def run_one(name):
    tree, model, ov, mode = RUNS[name]
    out = OUT / f"{name}.npz"
    js = out.with_suffix(".json")
    if js.exists():
        return name, json.loads(js.read_text())["status"]
    env = {"AUDIT_TAUMIN_MODE": mode} if mode is not None else {}
    st = harness.run_model(tree, P / f"{model}.ini", out, overrides=ov, threads=1, timeout=1500, extra_env=env)
    return name, st["status"]


def npz(name):
    return OUT / f"{name}.npz"


def summ(c, keys=("unlensed_scalar", "total")):
    return {k: {s: f"{v:.2e}" for s, v in c["cl"].get(k, {}).items()} for k in keys} | {
        "worst_pk": f"{c['worst_pk']:.2e}", "bit_identical": c["bit_identical"]}


def report():
    ref_ade = harness.AUDIT_ROOT / "runs" / "A0" / "5_ADE_1.npz"
    ref_quint = harness.AUDIT_ROOT / "runs" / "A0" / "5_quint_cubic.npz"
    st = {n: json.loads((OUT / f"{n}.json").read_text())["status"] for n in RUNS if (OUT / f"{n}.json").exists()}
    print("STATUS", json.dumps(st, indent=0))
    print("(a) bit-identity where A0 already runs")
    for n in ("fix0_ADE_default", "fix1_ADE_default", "fix2_ADE_default"):
        print(" ", n, "vs runs/A0/5_ADE_1:", summ(harness.compare(npz(n), ref_ade, lmax=3400)))
    for n in ("fix0_quint_default", "fix1_quint_default"):
        print(" ", n, "vs runs/A0/5_quint_cubic:", summ(harness.compare(npz(n), ref_quint, lmax=3400)))
    print("(b1) l_max_scalar=2000 (A0 refuses): mode1/mode2 vs A0 l_max 2500, ell<=1900; GR baseline same change")
    for n in ("fix1_ADE_l2000", "fix2_ADE_l2000"):
        if st.get(n) == "ok":
            print(" ", n, "vs A0_ADE_l2500:", summ(harness.compare(npz(n), npz("A0_ADE_l2500"), lmax=1900)))
    if st.get("fix1_ADE_l2000") == "ok" and st.get("fix2_ADE_l2000") == "ok":
        print("  mode1 vs mode2:", summ(harness.compare(npz("fix1_ADE_l2000"), npz("fix2_ADE_l2000"), lmax=1900)))
    print("  GR baseline A0_GR_l2000 vs A0_GR_l2500:",
          summ(harness.compare(npz("A0_GR_l2000"), npz("A0_GR_l2500"), lmax=1900)))
    print("(b2) no source windows, l_max 3500 (A0 refuses): mode1/2 vs A0 with windows, ell<=3000; GR baseline")
    for n in ("fix1_ADE_nowin", "fix2_ADE_nowin"):
        if st.get(n) == "ok":
            print(" ", n, "vs A0_ADE_win_notens:", summ(harness.compare(npz(n), npz("A0_ADE_win_notens"), lmax=3000)))
    if st.get("fix1_ADE_nowin") == "ok" and st.get("fix2_ADE_nowin") == "ok":
        print("  mode1 vs mode2:", summ(harness.compare(npz("fix1_ADE_nowin"), npz("fix2_ADE_nowin"), lmax=3000)))
    print("  GR baseline A0_GR_nowin vs A0_GR_win_notens:",
          summ(harness.compare(npz("A0_GR_nowin"), npz("A0_GR_win_notens"), lmax=3000)))
    print("(b3) transfer only, kmax=0.1 h/Mpc (A0 refuses): P(k) vs kmax=2 run on common k; GR baseline")
    for n in ("fix1_quint_tonly_k0.1", "fix2_quint_tonly_k0.1"):
        if st.get(n) == "ok":
            print(" ", n, "vs A0_quint_tonly_k2:", summ(harness.compare(npz(n), npz("A0_quint_tonly_k2")), ()))
    print("  GR baseline:", summ(harness.compare(npz("A0_GR_tonly_k0.1"), npz("A0_GR_tonly_k2")), ()))


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "run"
    if what == "run":
        names = sys.argv[2:] or list(RUNS)
        with ThreadPoolExecutor(4) as ex:
            for name, status in ex.map(run_one, names):
                print(name, status, flush=True)
    else:
        report()


if __name__ == "__main__":
    main()
