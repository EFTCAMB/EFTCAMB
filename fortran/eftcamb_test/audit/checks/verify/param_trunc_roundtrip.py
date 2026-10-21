"""Verifier (param-name-truncation): independent check through the Python wrapper (camb/eftcamb.py ->
11_EFTCAMB_python.f90 -> EFTCAMB_initialize -> model init_model_parameters_from_file).

Round trip: initialise a model with a given function-selection dict, ask EFTCAMB for its own parameter names
(EFTCAMB.param_names()), re-initialise passing every reported name with a distinct value, then read back
EFTCAMB.param_values().  Every name whose value does not come back is a key that EFTCAMB reports but never reads.

usage: python param_trunc_roundtrip.py [tree]      (tree = A0 by default, or a path to a tree copy)
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PY = sys.executable


def child(tree: str, cfg: dict) -> dict:
    os.environ["OMP_NUM_THREADS"] = "1"
    sys.path.insert(0, str(harness.tree_path(tree)))
    import camb

    def init(d):
        cp = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12)
        cp.EFTCAMB.initialize_parameters(cp, d, print_header=False)
        return cp  # keep the CAMBparams alive (cp.EFTCAMB is a view into its memory)

    base = dict(cfg)
    cp1 = init(base)
    names = cp1.EFTCAMB.param_names()
    vals = {n: round(0.1 + 0.01 * i, 4) for i, n in enumerate(names)}
    cp2 = init({**base, **vals})
    names2 = cp2.EFTCAMB.param_names()
    got = cp2.EFTCAMB.param_values()
    rows = []
    for n, g in zip(names2, got):
        rows.append({"name": n, "len": len(n), "set": vals.get(n), "got": g,
                     "ok": vals.get(n) is not None and abs(g - vals[n]) < 1e-12})
    return {"names_same": names == names2, "rows": rows}


CONFIGS = {}
_extra = {}
for fn in ("RPHmassP_ODE", "RPHkineticity_ODE", "RPHbraiding_ODE", "RPHtensor_ODE"):
    _extra.update({f"{fn}_Fourier_order_cos": 1, f"{fn}_Fourier_order_sin": 1, f"{fn}_Taylor_order": 2,
                   f"{fn}_Pade_order_N": 1, f"{fn}_Pade_order_D": 1})
for flag in range(1, 12):
    CONFIGS[f"RPH_flag{flag}"] = {"EFTflag": 2, "AltParEFTmodel": 1, "RPHwDE": 0, "RPHusealphaM": False,
                                  "RPHmassPmodel": 0, "RPHkineticitymodel": 0, "RPHbraidingmodel": 0,
                                  "RPHtensormodel": 0, "RPHmassPmodel_ODE": flag, "RPHkineticitymodel_ODE": flag,
                                  "RPHbraidingmodel_ODE": flag, "RPHtensormodel_ODE": flag, **_extra}
_extra3 = {}
for fn in ("Designer3p1sigma1_ODE", "Designer3p1rho1N_ODE", "Designer3p1rho2_ODE", "Designer3p1rho1_ODE"):
    _extra3.update({f"{fn}_Fourier_order_cos": 1, f"{fn}_Fourier_order_sin": 1, f"{fn}_Taylor_order": 2,
                    f"{fn}_Pade_order_N": 1, f"{fn}_Pade_order_D": 1})
for flag in range(1, 12):
    CONFIGS[f"3p1_flag{flag}"] = {"EFTflag": 3, "DesignerEFTmodel": 3, "EFTwDE": 0, "Designer3p1LuminalGW": False,
                                  "Designer3p1modelsigma1": 1, "Designer3p1modelrho1": 1,
                                  "Designer3p1modelrho1N": 0, "Designer3p1modelrho2": 0,
                                  "Designer3p1modelsigma1_ODE": flag, "Designer3p1modelrho1N_ODE": flag,
                                  "Designer3p1modelrho2_ODE": flag, "Designer3p1modelrho1_ODE": flag, **_extra3}


def main():
    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    out = {}
    for name, cfg in CONFIGS.items():
        with harness.slot(1):
            p = subprocess.run([PY, __file__, "--child", tree, json.dumps(cfg)], capture_output=True, text=True,
                               timeout=300, env=dict(os.environ, OMP_NUM_THREADS="1"))
        res = None
        for line in p.stdout.splitlines():
            if line.startswith("RESULT "):
                res = json.loads(line[7:])
        if res is None:
            out[name] = {"status": f"rc={p.returncode}", "tail": (p.stdout + p.stderr).splitlines()[-4:]}
            print(f"{name:14s} child failed rc={p.returncode}: {' / '.join(out[name]['tail'])}")
            continue
        out[name] = res
        bad = [r for r in res["rows"] if not r["ok"]]
        n_long = sum(r["len"] > 20 for r in res["rows"])
        print(f"{name:14s} params {len(res['rows']):3d}  names>20 chars {n_long:3d}  not read back {len(bad):3d}")
        for r in bad:
            print(f"      {r['name']:30s} len {r['len']:2d}  set {r['set']}  got {r['got']}")
        # sanity: every name <= 20 chars must round-trip, every name > 20 must fail
        wrong_short = [r["name"] for r in res["rows"] if r["len"] <= 20 and not r["ok"]]
        right_long = [r["name"] for r in res["rows"] if r["len"] > 20 and r["ok"]]
        if wrong_short or right_long:
            print(f"      !! not explained by length: short-not-read {wrong_short} long-read {right_long}")
    dest = Path(os.environ.get("VERIFY_OUT", harness.AUDIT_ROOT / "work" / "verify" / "param_trunc"))
    dest.mkdir(parents=True, exist_ok=True)
    tag = Path(str(tree)).name
    (dest / f"roundtrip_{tag}.json").write_text(json.dumps(out, indent=1))
    print("saved", dest / f"roundtrip_{tag}.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        r = child(sys.argv[2], json.loads(sys.argv[3]))
        print("RESULT " + json.dumps(r))
    else:
        main()
