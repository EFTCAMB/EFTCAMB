"""Independent check of finding get-eft-functions-no-cache-segfault.

Each case runs in its own child process (a crash only kills the child). The whole batch holds one machine-wide CPU
slot (harness.slot) and runs at nice 10.

Cases (tree selectable):
  gr_bg          GR (EFTflag=0) set_params, get_background, get_eft_functions
  gr_bg_nonu     same with mnu=0 (num_nu_massive=0): isolates the model dereference from the neutrino loop
  gr_full        GR, full get_results with Cls (lmax 500), get_eft_functions
  gr_evol        GR, full get_results, get_scale_evolution (EFTCAMB_GetEvolution path)
  eft_none_info  EFT ini, WantCls=WantTransfer=False, get_results: inspect res.Params (EFTflag/H0/allocated) only
  eft_none_evol  same, get_scale_evolution output (zeros claim)
  eft_none_fun   same, get_eft_functions
  eft_bg_fun     control: EFT ini, get_background, get_eft_functions
  eft_tr_fun     control: EFT ini, WantTransfer only, get_results, get_eft_functions
  b0_none_info   upstream B0: get_results with nothing to compute, inspect res.Params
"""

import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness  # noqa: E402

TREES = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees"
INI = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"
OUT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/geteft_noc"

PRE = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
case, ini = sys.argv[2], sys.argv[3]
a = np.array([0.1, 0.5, 1.0])
def gr(mnu=0.06, full=False):
    p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, mnu=mnu, lmax=500)
    p.WantCls = full
    p.WantTransfer = False
    return p
def eft():
    p = camb.read_ini(ini, no_validate=True)
    p.max_l = 500
    p.max_eta_k = 1500
    return p
def describe(res, label):
    P = res.Params
    try:
        ef = P.EFTCAMB
        efs = None if ef is None else int(ef.EFTflag)
    except Exception as e:
        efs = 'err:' + type(e).__name__ + ':' + str(e)[:80]
    try:
        ct = res.conformal_time(np.array([0.0, 1.0])).tolist()
    except Exception as e:
        ct = 'err:' + str(e)[:80]
    print('INFO', label, json.dumps({'H0': P.H0, 'ombh2': P.ombh2, 'EFTflag_in_results': efs, 'conformal_time(z=0,1)': ct}), flush=True)
"""

BODY = {
    "gr_bg": "p = gr(); res = camb.get_background(p); print('STEP before call', flush=True)\n"
             "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "gr_bg_nonu": "p = gr(mnu=0.0); print('num_nu_massive', p.num_nu_massive); res = camb.get_background(p)\n"
                  "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "gr_full": "p = gr(full=True); res = camb.get_results(p); print('STEP Cls done TT[100]=', res.get_cmb_power_spectra(CMB_unit='muK')['total'][100,0], flush=True)\n"
               "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "gr_evol": "p = gr(full=True); res = camb.get_results(p)\n"
               "f, r = p.EFTCAMB.get_scale_evolution(res, [0.01], a); print('RESULT evol adotoa', r['adotoa'].tolist(), flush=True)",
    "eft_none_info": "p = eft(); p.WantCls = False; p.WantTransfer = False; p.DoLensing = False\n"
                     "res = camb.get_results(p); print('params EFTflag', p.EFTCAMB.EFTflag); describe(res, 'A0-noout')",
    "eft_none_evol": "p = eft(); p.WantCls = False; p.WantTransfer = False; p.DoLensing = False\n"
                     "res = camb.get_results(p)\n"
                     "f, r = p.EFTCAMB.get_scale_evolution(res, [0.01, 0.1], a)\n"
                     "print('RESULT evol max|adotoa|', float(np.abs(r['adotoa']).max()), 'max|a|', float(np.abs(r['a']).max()), "
                     "'max|EFTOmegaV|', float(np.abs(r['EFTOmegaV']).max()), 'nonzero fields', "
                     "sum(1 for k in r.dtype.names if np.any(r[k] != 0)), 'of', len(r.dtype.names), flush=True)",
    "eft_none_fun": "p = eft(); p.WantCls = False; p.WantTransfer = False; p.DoLensing = False\n"
                    "res = camb.get_results(p)\n"
                    "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "eft_bg_fun": "p = eft(); res = camb.get_background(p); describe(res, 'bg')\n"
                  "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), 'EFTOmega', r['EFTOmegaV'].tolist(), flush=True)",
    "eft_tr_fun": "p = eft(); p.WantCls = False; p.WantTransfer = True; res = camb.get_results(p); print('STEP transfer done s8=', res.get_sigma8_0(), flush=True); describe(res, 'tr')\n"
                  "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), 'EFTOmega', r['EFTOmegaV'].tolist(), flush=True)",
    "eft_trcls_fun": "p = eft(); p.WantCls = True; p.WantTransfer = True; res = camb.get_results(p); print('STEP cls+transfer done', flush=True); describe(res, 'trcls')\n"
                     "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "eft_cls_fun": "p = eft(); p.WantCls = True; p.WantTransfer = False; res = camb.get_results(p); describe(res, 'cls')\n"
                   "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "eft_bgnothermo_fun": "p = eft(); res = camb.get_background(p, no_thermo=True); describe(res, 'bgnt')\n"
                   "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "gr_bgnothermo": "p = gr(); res = camb.get_background(p, no_thermo=True); print('STEP before call', flush=True)\n"
             "f, r = p.EFTCAMB.get_eft_functions(res, a); print('RESULT adotoa', r['adotoa'].tolist(), flush=True)",
    "b0_none_info": "p = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965)\n"
                    "p.WantCls = False; p.WantTransfer = False; p.DoLensing = False\n"
                    "res = camb.get_results(p); describe(res, 'B0-noout')",
}


def run_case(tree: str, case: str) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    code = PRE + "\n" + BODY[case]
    tpath = tree if tree.startswith("/") else f"{TREES}/{tree}"
    r = subprocess.run([sys.executable, "-c", code, tpath, case, INI], env=env, capture_output=True,
                       text=True, timeout=500, preexec_fn=lambda: os.nice(10))
    status = "ok" if r.returncode == 0 else (f"signal {-r.returncode}" if r.returncode < 0 else f"rc={r.returncode}")
    lines = [ln for ln in r.stdout.splitlines() if ln.startswith(("INFO", "RESULT", "STEP", "params", "num_nu"))]
    err = [ln for ln in r.stderr.splitlines() if ln.strip()][-25:]
    return {"tree": tree, "case": case, "status": status, "stdout": lines, "stderr_tail": err}


def main() -> None:
    plan = [tuple(x.rsplit(":", 1)) for x in sys.argv[1:]] if len(sys.argv) > 1 else [
        ("A0", "gr_bg"), ("A0", "gr_bg_nonu"), ("A0", "gr_full"), ("A0", "gr_evol"), ("A0", "eft_none_info"),
        ("A0", "eft_none_evol"), ("A0", "eft_none_fun"), ("A0", "eft_bg_fun"), ("A0", "eft_tr_fun"),
        ("B0", "b0_none_info"), ("B1", "gr_bg"), ("B1", "eft_none_info"), ("B1", "eft_none_fun"), ("B1", "eft_bg_fun"),
    ]
    os.makedirs(OUT, exist_ok=True)
    allres = []
    with harness.slot(1):
        for tree, case in plan:
            res = run_case(tree, case)
            allres.append(res)
            print(f"[{tree}:{case}] {res['status']}")
            for ln in res["stdout"]:
                print("    ", ln[:400])
            if res["status"] != "ok":
                for ln in res["stderr_tail"][-25:]:
                    print("     !", ln[:300])
    tag = "_".join(sorted({os.path.basename(t) for t, _ in plan}))
    with open(os.path.join(OUT, f"indep_{tag}.json"), "w") as f:
        json.dump(allres, f, indent=1)


if __name__ == "__main__":
    main()
