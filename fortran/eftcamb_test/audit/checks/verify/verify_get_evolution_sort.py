"""Independent adversarial check of finding get-evolution-sort-misalign.

Claim: EFTCAMB.get_evolution (camb/eftcamb.py) sorts the requested conformal times and returns rows in sorted order,
so for unsorted / decreasing input the rows do not match the caller's times.

Independent method (different from the repro):
  A. models different from the repro (PEFT, quintessence full mapping, mc5e designer), a random permutation of eta
     and an unsorted k list; compare the returned 'tau' and 'a' fields with
       - the requested eta (caller order) and sorted(eta),
       - an independent scale factor from results.redshift_at_conformal_time(eta) (CAMB background, not EFTCAMB),
       - upstream CAMB results.get_time_evolution(k, eta, ['a']) on the same results (sibling API, restores order).
  B. permutation invariance: output of permuted input == output of sorted input, bit for bit.
  C. k order is preserved (k field per row == requested k, unsorted k list).
  D. proposed fix implemented locally (monkeypatch, no repo edit): caller-order echo and bit-identical values vs the
     sorted call.

Usage: python verify_get_evolution_sort.py [--tree PATH]
"""

import argparse
import json
import os
import sys
import tempfile
from ctypes import POINTER, byref, c_double, c_int
from pathlib import Path

import numpy as np

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
OUT = Path(ROOT) / "work" / "verify"

ap = argparse.ArgumentParser()
ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
args = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, args.tree)
import camb  # noqa: E402
from camb import eftcamb as eftmod  # noqa: E402

print("camb from", camb.__file__)


def fixed_get_evolution(eft, camb_results, q, eta):
    """Proposed fix (argsort + inverse permutation), calling the same Fortran routine."""
    k = np.asfortranarray(np.atleast_1d(np.asarray(q, dtype=np.float64)))
    times = np.atleast_1d(np.asarray(eta, dtype=np.float64))
    order = np.argsort(times, kind="stable")
    times_s = np.asfortranarray(times[order])
    results = (eftmod.EFTCAMB_timestep_cache * times_s.shape[0] * k.shape[0])()
    eftmod.EFTCAMB_quantities_evolution(
        byref(camb_results),
        byref(c_int(k.shape[0])),
        k.ctypes.data_as(POINTER(c_double)),
        byref(c_int(times_s.shape[0])),
        times_s.ctypes.data_as(POINTER(c_double)),
        byref(results),
    )
    fields = results[0][0].variables()
    results = np.ctypeslib.as_array(results)
    unsorted = np.empty_like(results)
    unsorted[:, order] = results
    return fields, unsorted


def num_fields(arr):
    return [n for n in arr.dtype.names if np.issubdtype(arr.dtype[n], np.floating)]


checks = []


def check(name, value, thr, passed):
    checks.append({"name": name, "value": float(value), "threshold": thr, "status": "PASS" if passed else "FAIL"})
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {value:.3e} (threshold {thr})")


models = ["2_PEFT_Omega_lin_1", "5_quint_exponential", "4_Des_mc5e_wCDM_1", "5_hdsk_TG_1"]
rng = np.random.default_rng(12345)
tmp = tempfile.mkdtemp(prefix="vsort_", dir=str(OUT))
for model in models:
    print(f"\n=== {model}")
    ini = os.path.join(tmp, f"{model}.ini")
    Path(ini).write_text(
        f"DEFAULT({PARAMS / (model + '.ini')})\noutput_root = {tmp}/{model}\nEFTCAMB_write_background = F\n"
        "feedback_level = 0\n"
    )
    cp = camb.read_ini(ini, no_validate=True)
    res = camb.get_background(cp)
    z = np.array([0.0, 0.3, 1.0, 2.5, 6.0, 20.0, 100.0])
    eta_sorted_z = res.conformal_time(z)  # decreasing (z increasing)
    perm = rng.permutation(len(z))
    eta = eta_sorted_z[perm]
    kreq = np.array([0.2, 0.01, 0.05])
    _, ev = cp.EFTCAMB.get_evolution(res, kreq, eta)
    _, ev_s = cp.EFTCAMB.get_evolution(res, kreq, np.sort(eta))
    a_indep = 1.0 / (1.0 + res.redshift_at_conformal_time(eta))
    tau_ret = ev["tau"]
    a_ret = ev["a"]
    print("  requested eta   :", np.round(eta, 2).tolist())
    print("  returned tau[0] :", np.round(tau_ret[0], 2).tolist())
    print("  independent a   :", np.round(a_indep, 5).tolist())
    print("  returned a[0]   :", np.round(a_ret[0], 5).tolist())
    # A: caller order mismatch
    check(f"{model}: max|tau_ret - eta_req| [Mpc] (claim: >0)", np.max(np.abs(tau_ret - eta[None, :])), "==0 expected", np.max(np.abs(tau_ret - eta[None, :])) == 0)
    check(f"{model}: max|tau_ret - sort(eta_req)| [Mpc]", np.max(np.abs(tau_ret - np.sort(eta)[None, :])), "==0", np.max(np.abs(tau_ret - np.sort(eta)[None, :])) == 0)
    rel_a = np.max(np.abs(a_ret - a_indep[None, :]) / a_indep[None, :])
    check(f"{model}: max rel |a_ret - a_indep(eta_req)| (claim: >>1e-6)", rel_a, "<=1e-4", rel_a <= 1e-4)
    rel_a_s = np.max(np.abs(a_ret - np.sort(a_indep)[None, :]) / np.sort(a_indep)[None, :])
    check(f"{model}: max rel |a_ret - sort(a_indep)|", rel_a_s, "<=1e-4", rel_a_s <= 1e-4)
    # sibling upstream API
    try:
        camb_a = res.get_time_evolution(kreq, eta, ["a"])[:, :, 0]
        rel_c = np.max(np.abs(camb_a - a_indep[None, :]) / a_indep[None, :])
        check(f"{model}: upstream get_time_evolution 'a' vs a_indep(eta_req) (caller order)", rel_c, "<=1e-4", rel_c <= 1e-4)
        rel_ce = np.max(np.abs(a_ret - camb_a) / camb_a)
        check(f"{model}: EFTCAMB get_evolution 'a' vs upstream get_time_evolution 'a' (same eta)", rel_ce, "<=1e-4", rel_ce <= 1e-4)
    except Exception as exc:  # noqa: BLE001
        print("  upstream get_time_evolution failed:", repr(exc))
    # B: permutation invariance (outputs ignore caller order entirely)
    maxd = max(np.nanmax(np.abs(ev[f] - ev_s[f])) for f in num_fields(ev))
    check(f"{model}: permuted-input output vs sorted-input output (max abs diff all fields)", maxd, "==0 (claim)", maxd == 0)
    # C: k order kept
    kd = np.max(np.abs(ev["k"][:, 0] - kreq) / kreq) if "k" in ev.dtype.names else np.nan
    check(f"{model}: k field vs requested (unsorted) k, rel", kd, "<=1e-12", kd <= 1e-12)
    # D: proposed fix
    _, ev_f = fixed_get_evolution(cp.EFTCAMB, res, kreq, eta)
    check(f"{model}: FIX max|tau - eta_req|", np.max(np.abs(ev_f["tau"] - eta[None, :])), "==0", np.max(np.abs(ev_f["tau"] - eta[None, :])) == 0)
    rel_af = np.max(np.abs(ev_f["a"] - a_indep[None, :]) / a_indep[None, :])
    check(f"{model}: FIX max rel |a - a_indep(eta_req)|", rel_af, "<=1e-4", rel_af <= 1e-4)
    inv = np.argsort(perm)  # not used directly; compare fix with sorted call re-permuted
    order = np.argsort(eta, kind="stable")
    maxdf = max(np.nanmax(np.abs(ev_f[f][:, order] - ev_s[f])) for f in num_fields(ev))
    check(f"{model}: FIX values (reordered) vs current sorted-input call, max abs diff all fields", maxdf, "==0", maxdf == 0)
    _, ev_fs = fixed_get_evolution(cp.EFTCAMB, res, kreq, np.sort(eta))
    maxdfs = max(np.nanmax(np.abs(ev_fs[f] - ev_s[f])) for f in num_fields(ev))
    check(f"{model}: FIX on sorted input == current code on sorted input", maxdfs, "==0", maxdfs == 0)
    # E: get_scale_evolution decreasing a vs get_eft_functions (caller order) on adotoa
    a_dec = 1.0 / (1.0 + z)
    _, sev = cp.EFTCAMB.get_scale_evolution(res, [0.1], a_dec)
    _, efn = cp.EFTCAMB.get_eft_functions(res, a_dec)
    rel_h = np.max(np.abs(sev["adotoa"][0] - efn["adotoa"]) / np.abs(efn["adotoa"]))
    check(f"{model}: get_scale_evolution(a decreasing) adotoa vs get_eft_functions adotoa, rel (claim: large)", rel_h, "<=1e-3", rel_h <= 1e-3)
    _, sevf = fixed_get_evolution(cp.EFTCAMB, res, [0.1], res.conformal_time(1.0 / a_dec - 1.0))
    rel_hf = np.max(np.abs(sevf["adotoa"][0] - efn["adotoa"]) / np.abs(efn["adotoa"]))
    check(f"{model}: FIX scale evolution adotoa vs get_eft_functions adotoa, rel", rel_hf, "<=1e-3", rel_hf <= 1e-3)

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "verify_get_evolution_sort.json").write_text(json.dumps(checks, indent=1))
print("\nwritten", OUT / "verify_get_evolution_sort.json")
