"""Independent verification of finding wde-integral-unsupported-parametrizations.

Different code path from the t4a-param repro: no ini file, no harness worker.  Each case is a child Python process
that builds the parameters with camb.set_params(**dict) (EFTCAMB initialised from the dict through
EFTCAMB.initialize_parameters) and calls camb.get_background inside try/except CAMBError/finally.  A case "killed"
means the child printed neither a result nor an exception and never reached its finally block (Fortran ERROR STOP).

For the cases that run, H(z) and EFT functions (grhov_t, EFTOmegaV) are compared with the EFTwDE=1 reference that
describes the same constant w.

Models: designer f(R) (EFTflag=3, DesignerEFTmodel=1), shift-symmetric alpha_B (EFTflag=2, AltParEFTmodel=4) and,
as a control showing that the parametrizations themselves work, PureEFT (EFTflag=1) which integrates rho_DE
numerically instead of calling integral().

usage: python verify_wde_integral.py [tree]      (tree default A0)
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "wde_integral"

COSMO = {"H0": 70.0, "ombh2": 0.0226, "omch2": 0.112, "mnu": 0.06, "tau": 0.09, "feedback_level": 1}
STAB = {"EFT_ghost_math_stability": False, "EFT_mass_math_stability": False, "EFT_ghost_stability": True,
        "EFT_gradient_stability": True, "EFT_mass_stability": False, "EFT_mass_stability_rate": 10.0,
        "EFT_additional_priors": True, "EFTCAMB_turn_on_time": 0.01, "EFTCAMB_stability_time": 1.0e-10,
        "EFTCAMB_stability_threshold": 0.0}


def wde_variants(w, extra_nonconst=True):
    """EFTwDE options all describing constant w (plus a few genuinely non-constant ones)."""
    v = {
        "wDE=1 const": {"EFTwDE": 1, "EFTw0": w},
        "wDE=8 steplog v1=v2": {"EFTwDE": 8, "EFTw_v1": w, "EFTw_v2": w, "EFTw_at": 0.5, "EFTw_delta": 0.3},
        "wDE=6 Pade N=D=0": {"EFTwDE": 6, "EFTw_Pade_order_N": 0, "EFTw_Pade_order_D": 0, "EFTwa0": 0.0,
                             "EFTwN0": w},
        "wDE=6 Pade N=D=1 (const)": {"EFTwDE": 6, "EFTw_Pade_order_N": 1, "EFTw_Pade_order_D": 1, "EFTwa0": 1.0,
                                     "EFTwN0": w, "EFTwN1": 0.5 * w, "EFTwD1": 0.5},
        "wDE=7 Fourier order 0": {"EFTwDE": 7, "EFTw_a0": w, "EFTw_phase": 0.0},
        "wDE=11 exp2 lambda=0": {"EFTwDE": 11, "EFTw0": w, "EFTw1": 0.0},
    }
    if extra_nonconst:
        v["wDE=11 exp2 lambda=0.2 (non-const)"] = {"EFTwDE": 11, "EFTw0": w, "EFTw1": 0.2}
        v["wDE=7 Fourier cos1 (non-const)"] = {"EFTwDE": 7, "EFTw_Fourier_order_cos": 1, "EFTw_a0": w,
                                               "EFTw_phase": 0.0, "EFTw_a1": 0.02}
    return v


MODELS = {
    "DesfR": ({"EFTflag": 3, "DesignerEFTmodel": 1, "EFTB0": 1.0}, -0.9),
    "ShiftSym": ({"EFTflag": 2, "AltParEFTmodel": 4, "Shift_Symmetric_alphaB0": 0.4, "Shift_Symmetric_m": 1.0,
                  "Shift_Symmetric_alphaK0": 1.0}, -0.95),
    "PureEFT": ({"EFTflag": 1, "PureEFTmodel": 1, "PureEFTHorndeski": False,
                 **{f"PureEFTmodel{f}": 0 for f in ("Omega", "Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5",
                                                    "Gamma6")}}, -0.8),
}

ZS = [0.0, 0.5, 1.0, 2.0, 5.0]
AS = [10.0 ** (-3 + 3 * i / 19) for i in range(20)]


def child(tree, params):
    sys.path.insert(0, str(tree))
    import numpy as np

    import camb
    from camb.baseconfig import CAMBError

    print("STAGE start", flush=True)
    try:
        pars = camb.set_params(**params)
        print("STAGE set_params done", flush=True)
        res = camb.get_background(pars)
        print("STAGE get_background done", flush=True)
        H = [float(np.ravel(res.hubble_parameter(z))[0]) for z in ZS]
        fields, arr = pars.EFTCAMB.get_eft_functions(res, np.array(AS))
        out = {"H": H}
        for f in ("grhov_t", "EFTOmegaV", "adotoa"):
            out[f] = [float(x) for x in arr[f]]
        print("RESULT " + json.dumps(out), flush=True)
    except CAMBError as e:
        print("CAUGHT CAMBError: " + " | ".join(str(e).splitlines()[:3]), flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"CAUGHT {type(e).__name__}: {e}", flush=True)
    finally:
        print("FINALLY_REACHED", flush=True)


def run_case(tree, params, tag):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    cmd = [sys.executable, str(Path(__file__).resolve()), "--child", str(tree), json.dumps(params)]
    with harness.slot(1):
        p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600, cwd=OUT)
    log = OUT / f"{tag}.log"
    log.write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    lines = p.stdout.splitlines()
    res = next((json.loads(x[7:]) for x in lines if x.startswith("RESULT ")), None)
    caught = next((x for x in lines if x.startswith("CAUGHT")), None)
    fin = "FINALLY_REACHED" in lines
    err_stop = "ERROR STOP" in (p.stdout + p.stderr)
    notimpl = next((x.strip() for x in lines if "not implemented" in x), None)
    return {"rc": p.returncode, "result": res, "caught": caught, "finally": fin, "error_stop": err_stop,
            "not_implemented_msg": notimpl, "log": str(log)}


def maxrel(a, b):
    return max(abs(x - y) / max(abs(y), 1e-300) for x, y in zip(a, b))


def main():
    tree = harness.tree_path(sys.argv[1] if len(sys.argv) > 1 else "A0")
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {}
    for mname, (base, w) in MODELS.items():
        ref = None
        for iv, (vname, wd) in enumerate(wde_variants(w, extra_nonconst=(mname != "PureEFT")).items()):
            params = {**COSMO, **STAB, **base, **wd}
            tag = f"{mname}_{iv}_{vname.split()[0].replace('=', '')}"
            r = run_case(tree, params, tag)
            if r["result"] is not None:
                if ref is None and vname.startswith("wDE=1"):
                    ref = r["result"]
                    outcome = "ok (reference)"
                elif "non-const" in vname:
                    outcome = "ok (non-constant w, not compared)"
                else:
                    d = {k: maxrel(r["result"][k], ref[k]) for k in ("H", "grhov_t", "EFTOmegaV")}
                    r["maxrel_vs_ref"] = d
                    outcome = "ok, max rel diff vs ref: " + ", ".join(f"{k}={v:.2e}" for k, v in d.items())
            elif r["caught"]:
                outcome = f"python exception caught (rc={r['rc']}): {r['caught'][:150]}"
            else:
                outcome = (f"KILLED rc={r['rc']} error_stop={r['error_stop']} finally_reached={r['finally']} "
                           f"msg='{r['not_implemented_msg']}'")
            print(f"{mname:9s} w={w:5.2f} {vname:38s} -> {outcome}", flush=True)
            summary[f"{mname}|{vname}"] = r | {"outcome": outcome}
    (OUT / f"summary_{Path(tree).name}.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        child(sys.argv[2], json.loads(sys.argv[3]))
    else:
        main()
