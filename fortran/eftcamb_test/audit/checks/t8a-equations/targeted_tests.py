"""Targeted flag tests on the audit-patched tree P (see make_patched_tree.py).

  polter  : EFTflag=0, omk=-0.05/+0.05: P with AUDIT_POLTER_GR=1 vs plain CAMB (B0 runs from gr_consistency.py)
  dpmatter: 5_hdsk_TG_* (evolve_delta_phi) with do_late_rad_truncation=T, omnuh2=0.01: flag off vs on
  metric_h: 5_hdsk_TG_1 with EFTCAMB_evolve_metric_h=T: flag off vs AUDIT_FIX_METRIC_H
  dgp     : 5_hdsk_TG_1/3 (evolve_delta_phi) with omnuh2=0.01: massive-nu dp = p_nu*delta_nu (code) vs exact dp
            (AUDIT_DGP_EXACT=1); reuses the dpmatter "off" runs when present
Usage: python targeted_tests.py WHICH
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8a-equations")
P = W / "P"
OUT = W / "targeted"
OUT.mkdir(exist_ok=True)


def summary(c):
    return {"bit": c["bit_identical"], "worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"],
            "unlensed_scalar": c["cl"].get("unlensed_scalar"), "tensor": c["cl"].get("tensor"),
            "sigma8": c.get("sigma8")}


def polter():
    res = {}
    base = harness.PARAMS_DIR / "base_params.ini"
    for name, omk in (("omk-0.05", "-0.05"), ("omk+0.05", "0.05"), ("flat", "0")):
        ov = {"EFTflag": "0", "get_tensor_cls": "F", "num_redshiftwindows": "0", "omk": omk}
        if omk == "0":
            ov.pop("omk")
        r = harness.run_model(P, base, OUT / f"polter_{name}.npz", overrides=ov,
                              python_set={"Reion.max_redshift": 50.0}, extra_env={"AUDIT_POLTER_GR": "1"})
        b = W / "gr" / f"B0_{name}.npz"
        res[name] = summary(harness.compare(OUT / f"polter_{name}.npz", b)) if r["status"] == "ok" else r["status"]
    return res


def dpmatter():
    res = {}
    ov = {"do_late_rad_truncation": "T", "omnuh2": "0.01"}
    jobs = []
    for m in ("5_hdsk_TG_1", "5_hdsk_TG_3"):
        ini = harness.PARAMS_DIR / f"{m}.ini"
        jobs.append((m, "off", ini, {}))
        jobs.append((m, "on", ini, {"AUDIT_FIX_DPMATTER": "1"}))
    with ThreadPoolExecutor(4) as ex:
        st = list(ex.map(lambda j: harness.run_model(P, j[2], OUT / f"dpm_{j[0]}_{j[1]}.npz", overrides=ov,
                                                     extra_env=j[3]), jobs))
    for m in ("5_hdsk_TG_1", "5_hdsk_TG_3"):
        try:
            res[m] = summary(harness.compare(OUT / f"dpm_{m}_on.npz", OUT / f"dpm_{m}_off.npz"))
        except Exception as e:  # noqa: BLE001
            res[m] = repr(e) + str([s["status"] for s in st])
    return res


def metric_h():
    res = {}
    ini = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
    ov = {"EFTCAMB_evolve_metric_h": "T"}
    jobs = [("off", {}), ("on", {"AUDIT_FIX_METRIC_H": "1"})]
    with ThreadPoolExecutor(2) as ex:
        st = list(ex.map(lambda j: harness.run_model(P, ini, OUT / f"mh_{j[0]}.npz", overrides=ov, extra_env=j[1]),
                         jobs))
    res["status"] = [s["status"] for s in st]
    if all(s["status"] == "ok" for s in st):
        res["on_vs_off"] = summary(harness.compare(OUT / "mh_on.npz", OUT / "mh_off.npz"))
        ref = harness.PARAMS_DIR.parent / "audit"  # placeholder to keep paths explicit
        a0 = Path(str(harness.AUDIT_ROOT / "runs" / "A0" / "5_hdsk_TG_1.npz"))
        res["metric_h_vs_constraint_A0"] = summary(harness.compare(OUT / "mh_off.npz", a0))
        del ref
    return res


def dgp():
    res = {}
    ov = {"do_late_rad_truncation": "T", "omnuh2": "0.01"}
    jobs = []
    for m in ("5_hdsk_TG_1", "5_hdsk_TG_3"):
        ini = harness.PARAMS_DIR / f"{m}.ini"
        if not (OUT / f"dpm_{m}_off.npz").exists():
            jobs.append((ini, OUT / f"dpm_{m}_off.npz", {}))
        jobs.append((ini, OUT / f"dgp_{m}_exact.npz", {"AUDIT_DGP_EXACT": "1"}))
    with ThreadPoolExecutor(4) as ex:
        st = list(ex.map(lambda j: harness.run_model(P, j[0], j[1], overrides=ov, extra_env=j[2]), jobs))
    for m in ("5_hdsk_TG_1", "5_hdsk_TG_3"):
        try:
            res[m] = summary(harness.compare(OUT / f"dgp_{m}_exact.npz", OUT / f"dpm_{m}_off.npz"))
        except Exception as e:  # noqa: BLE001
            res[m] = repr(e) + str([s["status"] for s in st])
    return res


if __name__ == "__main__":
    which = sys.argv[1]
    r = {"polter": polter, "dpmatter": dpmatter, "metric_h": metric_h, "dgp": dgp}[which]()
    print(json.dumps(r, indent=1))
    (OUT / f"{which}.json").write_text(json.dumps(r, indent=1))
