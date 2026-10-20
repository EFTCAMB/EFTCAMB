"""Verifier (t3-eft-turnon-before-tauminn-lmax): run early-turn-on models in ordinary configurations
(no extra source windows; CMB-only at various lmax; transfer-only) on a tree and report pass/refusal and the
tauminn quoted in the refusal, next to the prediction tauminn = 0.95*min(0.001/maxq, 0.1),
maxq = max(2.5*lmax/tau0 [CMB], 1.05*kmax [transfer]).

Usage: python tauminn_configs.py <tree> <tag> [case ...]   (default: all cases)
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "tauminn"
TAU0 = 13991.2  # Mpc, from results.tau0 for these base params (A0)
H = 0.70
NOWIN = {"num_redshiftwindows": "0"}
NOTENS = {"get_tensor_cls": "F"}
TONLY = {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_vector_cls": "F", "do_lensing": "F", "num_redshiftwindows": "0"}

CASES = {
    # name: (model, overrides, lmax or None, transfer kmax in h/Mpc or None)
    "ADE_default": ("5_ADE_1", {}, None, None),
    "ADE_nowin_l3500": ("5_ADE_1", {**NOWIN, **NOTENS}, 3500, 2.0),
    "ADE_nowin_l2500": ("5_ADE_1", {**NOWIN, **NOTENS, "l_max_scalar": "2500"}, 2500, 2.0),
    "ADE_tonly_k2": ("5_ADE_1", TONLY, None, 2.0),
    "quintcubic_nowin_l1000": ("5_quint_cubic", {**NOWIN, **NOTENS, "l_max_scalar": "1000"}, 1000, 2.0),
    "quintcubic_nowin_l1500": ("5_quint_cubic", {**NOWIN, **NOTENS, "l_max_scalar": "1500"}, 1500, 2.0),
    "quintcubic_tonly_k0.1": ("5_quint_cubic", {**TONLY, "transfer_kmax": "0.1"}, None, 0.1),
    "quintcubic_tonly_k2": ("5_quint_cubic", TONLY, None, 2.0),
    "quintcubic_default": ("5_quint_cubic", {}, None, None),
    "GR_nowin_l1000": ("base_params", {**NOWIN, **NOTENS, "l_max_scalar": "1000"}, 1000, 2.0),
}


def predicted_tauminn(lmax, kmax_h):
    maxq = 0.0
    if lmax:
        maxq = 2.5 * lmax / TAU0
    if kmax_h:
        maxq = max(maxq, 1.05 * kmax_h * H)
    return 0.95 * min(0.001 / maxq, 0.1) if maxq else None


def job(tree, tag, name):
    model, ov, lmax, kmax = CASES[name]
    ini = harness.PARAMS_DIR / f"{model}.ini"
    out = OUT / f"{tag}_{name}.npz"
    st = harness.run_model(tree, ini, out, overrides=ov, threads=1, timeout=1500)
    log = out.with_suffix(".log").read_text(errors="replace")
    msg = re.findall(r"EFTCAMB starts before thermo tauminn[^\n]*", log)
    return {"case": name, "status": st["status"], "refusal": msg[:1],
            "predicted_tauminn": predicted_tauminn(lmax, kmax)}


def main():
    tree, tag = sys.argv[1], sys.argv[2]
    names = sys.argv[3:] or list(CASES)
    with ThreadPoolExecutor(4) as ex:
        for r in ex.map(lambda n: job(tree, tag, n), names):
            print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main()
