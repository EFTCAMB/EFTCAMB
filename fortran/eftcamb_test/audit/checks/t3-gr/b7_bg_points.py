"""t3-gr (b): K-mouflage / K-mimic / ADE do not converge to LCDM as the coupling -> 0 (H floor ~3-7e-5, unlensed EE
~1-2e-4). Hypothesis: the background tables (model_background_num_points = 1000 / 1500, equispaced in ln a, linear
interpolation) are too coarse. Rerun the near-GR points with 10x more background points.
"""

from concurrent.futures import ThreadPoolExecutor

from common import PARAMS, RUNS, harness

SUB = RUNS / "b7"
B1 = RUNS / "b1"
JOBS = {
    "Kmouflage_s1em08_np10000": (PARAMS / "5_Kmouflage_1.ini", {"eps2_0": "-1e-08", "model_background_num_points": "10000"},
                                 "REF_nu0__A0g", "Kmouflage__s1em08"),
    "Kmimic_s1em08_np10000": (PARAMS / "5_Kmimic_1.ini", {"eps2_0": "1e-08", "model_background_num_points": "10000"},
                              "REF_nu0__A0g", "Kmimic__s1em08"),
    "ADE_s1em08_np15000": (PARAMS / "5_ADE_1.ini", {"f_ac": "1e-08", "model_background_num_points": "15000"},
                           "REF_base__A0g", "ADE__s1em08"),
}


def one(item):
    name, (ini, ov, ref, base) = item
    out = SUB / f"{name}.npz"
    if not out.exists():
        harness.run_model("A0", ini, out, overrides=ov, threads=1, timeout=1500)
    return name, out, ref, base


def brief(p, ref):
    r = harness.compare(p, B1 / f"{ref}.npz")
    u = r["cl"]["unlensed_scalar"]
    return (f"uTT={u['TT']:.1e} uEE={u['EE']:.1e} pk={r['worst_pk']:.1e} s8={r['sigma8']:.1e} H={r['bg_H']:.1e}")


with ThreadPoolExecutor(3) as ex:
    for name, out, ref, base in ex.map(one, JOBS.items()):
        if not out.exists():
            print(name, "FAILED")
            continue
        print(f"{base:22s} default points: {brief(B1 / f'{base}.npz', ref)}")
        print(f"{name:22s}               : {brief(out, ref)}")
