"""Track 9 (v): error handling of the python wrapper.

Each case runs in its own subprocess (tree A0):
  1. performs an invalid action (bad model selection number, missing parameter, unstable model, NaN/inf/negative
     parameters, ...) inside try/except and records the Python exception (type + message), or "no exception";
  2. if the process survived: runs a reference model (2_PEFT_Omega_power_1 with reduced lmax) in the same process
     and saves its Cls / sigma8; these are compared to the same reference computed in a fresh process.
Variants "<case>+reuse" call get_results on the parameters object left behind by the failed call.

Usage: python probe_errors.py [case ...]
"""

from __future__ import annotations

import json
import sys

import numpy as np

from t9lib import WORK, run_snippet

OUT = WORK / "errors"
REF_INI = OUT / "ref_fast.ini"

REF = r"""
def reference(tag):
    p = camb.read_ini(str(WORK / 'errors' / 'ref_fast.ini'))
    r = camb.get_results(p)
    cl = r.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)
    np.savez(str(WORK / 'errors' / (tag + '.npz')), total=cl['total'], lens=cl['lens_potential'],
             sigma8=np.atleast_1d(r.get_sigma8()))
    emit(reference_done=True)

def base_dict(**kw):
    d = {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1, 'EFTOmega0': 0.01, 'EFTwDE': 0,
         'PureEFTmodelGamma1': 0, 'PureEFTmodelGamma2': 0, 'PureEFTmodelGamma3': 0, 'PureEFTmodelGamma4': 0,
         'PureEFTmodelGamma5': 0, 'PureEFTmodelGamma6': 0}
    d.update(kw)
    return d

def attempt(fn):
    try:
        v = fn()
        emit(exception=None, returned=repr(v)[:200])
        return True
    except Exception as e:
        emit(exception=type(e).__name__, message=str(e)[:300])
        return False

def fresh_pars():
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    p.set_for_lmax(1000, lens_potential_accuracy=0)
    return p
"""

PARAMS_DIR = "PARAMS"

CASES: dict[str, str] = {
    "fresh_reference": "reference('fresh')",
    # invalid model selection numbers through initialize_parameters
    "bad_EFTflag": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, base_dict(EFTflag=9), "
    "print_header=False))",
    "bad_PureEFTmodel": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, "
    "base_dict(PureEFTmodel=9), print_header=False))",
    "bad_Omega_param": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, "
    "base_dict(PureEFTmodelOmega=99), print_header=False))",
    "bad_EFTwDE": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, base_dict(EFTwDE=99), "
    "print_header=False))",
    "bad_Designer": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, {'EFTflag': 3, "
    "'DesignerEFTmodel': 9}, print_header=False))",
    "bad_FullMapping": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, {'EFTflag': 4, "
    "'FullMappingEFTmodel': 99}, print_header=False))",
    # missing parameter: constant Omega without EFTOmega0
    "missing_EFTOmega0": "p = fresh_pars(); d = base_dict(); d.pop('EFTOmega0'); "
    "attempt(lambda: p.EFTCAMB.initialize_parameters(p, d, print_header=False)); "
    "emit(read=p.EFTCAMB.read_parameters().get('EFTOmega0'), values=p.EFTCAMB.param_values())",
    # wrong types
    "string_value": "p = fresh_pars(); attempt(lambda: p.EFTCAMB.initialize_parameters(p, base_dict(EFTOmega0='abc'), "
    "print_header=False))",
    # read_ini with invalid selection (goes through validate_ini_file)
    "read_ini_bad_EFTflag": "f = WORK / 'errors' / 'bad_flag.ini'; "
    "f.write_text(f'DEFAULT({PARAMS}/2_PEFT_Omega_power_1.ini)\\nEFTflag = 9\\n'); "
    "attempt(lambda: camb.read_ini(str(f)))",
    "read_ini_bad_param": "f = WORK / 'errors' / 'bad_param.ini'; "
    "f.write_text(f'DEFAULT({PARAMS}/2_PEFT_Omega_power_1.ini)\\nPureEFTmodelOmega = 99\\n'); "
    "attempt(lambda: camb.read_ini(str(f)))",
    # unstable model (ghost): Omega const = -0.5
    "unstable": "p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(EFTOmega0=-0.5), print_header=False); "
    "attempt(lambda: camb.get_results(p))",
    "unstable_background": "p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(EFTOmega0=-0.5), "
    "print_header=False); attempt(lambda: camb.get_background(p))",
    # NaN / inf / negative
    "nan_param": "p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(EFTOmega0=float('nan')), "
    "print_header=False); emit(values=p.EFTCAMB.param_values()); "
    "ok = attempt(lambda: camb.get_results(p)); "
    "r = camb.get_results(p) if ok else None; "
    "emit(cl_finite=bool(np.isfinite(r.get_cmb_power_spectra(raw_cl=True)['total']).all()) if r else None)",
    "inf_param": "p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(EFTOmega0=float('inf')), "
    "print_header=False); emit(values=p.EFTCAMB.param_values()); attempt(lambda: camb.get_results(p))",
    "nan_turn_on": "p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(EFTCAMB_turn_on_time=float('nan')),"
    " print_header=False); attempt(lambda: (camb.get_results(p), None)[1])",
    "negative_turn_on": "p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(EFTCAMB_turn_on_time=-0.1),"
    " print_header=False); attempt(lambda: (camb.get_results(p), None)[1])",
    "negative_bg_points": f"p = camb.read_ini(str({PARAMS_DIR} / '4_Des_F_R_LCDM_1.ini')); "
    "d = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in p.EFTCAMB.read_parameters().items()}; "
    "d['model_background_num_points'] = -10; "
    "attempt(lambda: p.EFTCAMB.initialize_parameters(p, d, print_header=False)); attempt(lambda: camb.get_results(p))",
}
# controls: a valid run of another model / of a GR model / of the same reference before the reference
CASES["control_other_model_first"] = (f"q = camb.read_ini(str({PARAMS_DIR} / '4_Des_F_R_LCDM_1.ini')); "
                                      "attempt(lambda: (camb.get_results(q), None)[1])")
CASES["control_gr_lmax2500_first"] = ("q = fresh_pars(); q.set_for_lmax(2500, lens_potential_accuracy=1); "
                                      "attempt(lambda: (camb.get_results(q), None)[1])")
CASES["control_reference_twice"] = "reference('control_reference_twice_first')"
CASES["control_fresh_eft_first"] = ("p = fresh_pars(); p.EFTCAMB.initialize_parameters(p, base_dict(), "
                                  "print_header=False); attempt(lambda: (camb.get_results(p), None)[1])")
# re-use of the parameters left behind by a failed initialization
for c in ("bad_EFTflag", "bad_PureEFTmodel", "bad_Omega_param", "bad_Designer"):
    CASES[c + "+reuse"] = CASES[c] + "; emit(EFTflag_left=p.EFTCAMB.EFTflag); attempt(lambda: camb.get_results(p))"
for c in list(CASES):
    if c != "fresh_reference" and "+reuse" not in c:
        CASES[c] = CASES[c] + f"; reference({c!r})"


def compare(tag: str) -> dict:
    a = np.load(OUT / "fresh.npz")
    f = OUT / f"{tag}.npz"
    if not f.exists():
        return {"after_reference": "missing"}
    b = np.load(f)
    same = all(np.array_equal(a[k], b[k]) for k in a.files)
    rel = max(float(np.max(np.abs(a[k] - b[k]) / np.maximum(np.abs(a[k]), 1e-300))) for k in a.files)
    return {"after_identical_to_fresh": same, "after_maxrel": rel}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    from t9lib import PARAMS

    REF_INI.write_text(f"DEFAULT({PARAMS}/2_PEFT_Omega_power_1.ini)\nl_max_scalar = 1200\nget_tensor_cls = F\n"
                       "EFTCAMB_write_background = F\n")
    sel = sys.argv[1:] or list(CASES)
    if "fresh_reference" not in sel and not (OUT / "fresh.npz").exists():
        sel = ["fresh_reference"] + sel
    summary = {}
    for name in sel:
        r = run_snippet(REF + CASES[name], threads=2, timeout=900, cwd=str(OUT))
        rec = {"status": r["status"], **r["result"]}
        if rec.get("reference_done") and name != "fresh_reference":
            rec.update(compare(name))
        if r["status"] != "ok" or rec.get("exception") is None:
            rec["tail"] = [t for t in r["tail"] if t.strip()][-4:]
        summary[name] = rec
        print(f"== {name}: {json.dumps(rec, default=str)}", flush=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
