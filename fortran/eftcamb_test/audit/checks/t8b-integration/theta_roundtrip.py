"""t8b check: H0 <-> cosmomc_theta round trip through camb/model.py set_cosmology (EFTCAMB MOD min_H0 bracket)
and the EFT-routed CosmomcTheta, for every test model. Background only.
For each model: set_cosmology(H0=H0_ini, ...) -> theta = cosmomc_theta(); set_cosmology(cosmomc_theta=theta, ...)
-> H0_rec. Expected |H0_rec/H0_ini - 1| < 1e-4 (brentq rtol 5e-5). Usage: python theta_roundtrip.py TREE [pattern]"""

import sys
import time


def main(args):
    import camb

    from t8b_common import read_model

    out = {}
    for ini in args["inis"]:
        t0 = time.time()
        try:
            cp = read_model(camb, ini)
            kw = dict(ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, omnuh2_active=cp.omnuh2, mnu=None,
                      num_massive_neutrinos=max(1, int(cp.num_nu_massive)) if cp.omnuh2 > 0 else 0,
                      nnu=cp.N_eff, TCMB=cp.TCMB, YHe=cp.YHe)
            if cp.omnuh2 == 0:
                kw["omnuh2_active"] = 0.0
            H0 = cp.H0
            cp.set_cosmology(H0=H0, **kw)
            res = camb.get_background(cp, no_thermo=True)
            th = res.cosmomc_theta()
            cp.set_cosmology(cosmomc_theta=th, **kw)
            out[ini] = {"H0_in": H0, "theta": th, "H0_rec": cp.H0, "rel": cp.H0 / H0 - 1,
                        "t": round(time.time() - t0, 2)}
        except Exception as e:  # noqa: BLE001
            out[ini] = {"error": repr(e)[:200], "t": round(time.time() - t0, 2)}
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, PARAMS, WORK, run_in_tree

    tree = sys.argv[1]
    pattern = sys.argv[2] if len(sys.argv) > 2 else "*.ini"
    inis = sorted(p.name for p in PARAMS.glob(pattern) if "base_params" not in p.name)
    res = {}
    chunk = 8
    for i in range(0, len(inis), chunk):
        r = run_in_tree(tree, HERE / "theta_roundtrip.py", {"inis": inis[i:i + chunk]}, timeout=1200)
        if "error" in r and "stdout_tail" in r:
            # crashed: run one by one
            for ini in inis[i:i + chunk]:
                r1 = run_in_tree(tree, HERE / "theta_roundtrip.py", {"inis": [ini]}, timeout=600)
                res.update(r1 if ini in r1 else {ini: {"error": "crash " + r1.get("error", "")}})
        else:
            res.update(r)
    name = str(tree).rstrip("/").split("/")[-1]
    (WORK / "out" / f"theta_roundtrip_{name}.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        if "error" in v:
            print(f"{k:32s} ERROR {v['error']}")
        else:
            print(f"{k:32s} H0_in={v['H0_in']:.4f} H0_rec={v['H0_rec']:.4f} rel={v['rel']:+.2e} t={v['t']}")
