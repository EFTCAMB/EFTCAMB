"""CAMB-level impact of finding spline5-poly5-dx6: A0 vs a tree copy with ONLY the one-line poly5 fix.

Model: pure EFT, Omega = 1.0 a^2 (2_PEFT_Omega_power_1.ini), Omega written as a quintic spline (PureEFTmodelOmega=10):
  logx_wide : 9 nodes a = 1e-4 ... 1 (half-decade steps), Logx=T  (interval in ln a = 1.15 > 1)
  logx_mod  : 8 nodes a = 0.01 ... 1, Logx=T
  lin_ctrl  : 6 nodes a = 0 ... 1 linear (a^2 is a quadratic, so c6 = 0: negative control, no change expected)
Reference: power-law run of the same model (runs/A0/2_PEFT_Omega_power_1.npz).

usage: python verify_spline5_cl_impact.py build|run|compare
"""

import json
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "spline5_poly5"
FIXTREE = WORK / "fix_poly5"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_power_1.ini"
REF = harness.AUDIT_ROOT / "runs" / "A0" / "2_PEFT_Omega_power_1.npz"
OLD = "    poly5 = c(1) + c(2)*dx + c(3)*dx**2 + c(4)*dx**3 + c(5)*dx**4 + c(6)*dx**6"
NEW = "    poly5 = c(1) + c(2)*dx + c(3)*dx**2 + c(4)*dx**3 + c(5)*dx**4 + c(6)*dx**5"


def spl(xs, vs, logx):
    d = {"PureEFTmodelOmega": 10, "EFTOmega_Spline_Order": 5, "EFTOmega_Logx_Flag": "T" if logx else "F",
         "EFTOmega_Spline_Pixels": len(xs)}
    for i, (x, v) in enumerate(zip(xs, vs)):
        d[f"EFTOmegax{i + 1}"] = repr(float(x))
        d[f"EFTOmegav{i + 1}"] = repr(float(v))
    return d


WIDE = [10.0 ** (-4 + 0.5 * i) for i in range(9)]
MOD = [0.01, 0.03, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0]
LIN = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
LIN11 = [0.1 * i for i in range(11)]
LOG15 = [10.0 ** (-3 + 3 * i / 14) for i in range(15)]
COARSE = [0.01, 0.05, 0.2, 0.5, 0.75, 1.0]
CASES = {
    # first attempt (both A0 and FIX: ghost instability near the first node -> not usable, kept for the record)
    "logx_wide": spl(WIDE, [a**2 for a in WIDE], True),
    "logx_mod": spl(MOD, [a**2 for a in MOD], True),
    "lin_ctrl": spl(LIN, [a**2 for a in LIN], False),
    # second attempt
    "lin_pow25": spl(LIN, [a**2.5 for a in LIN], False),
    "lin11_pow25": spl(LIN11, [a**2.5 for a in LIN11], False),
    "log15_a2": spl(LOG15, [a**2 for a in LOG15], True),
    "coarse_a2": spl(COARSE, [a**2 for a in COARSE], True),
}
REFS = {"lin_pow25": "pow25_ref"}


def build():
    if not FIXTREE.exists():
        shutil.copytree(harness.TREES["A0"], FIXTREE, symlinks=True)
    f = FIXTREE / "fortran/eftcamb/02_bundled_algorithms/02p12_splines5.f90"
    s = f.read_text()
    if OLD in s:
        f.write_text(s.replace(OLD, NEW))
    assert NEW in f.read_text()
    print(harness.build_tree("A0", FIXTREE))


def run(only=""):
    jobs = [(t, n) for t in ("A0", "FIX") for n in CASES if (not only or n in only.split(","))]

    def one(job):
        t, n = job
        tree = "A0" if t == "A0" else FIXTREE
        return job, harness.run_model(tree, INI, WORK / t / f"{n}.npz", overrides=CASES[n], threads=1, timeout=1500)

    with ThreadPoolExecutor(4) as ex:
        for job, st in ex.map(one, jobs):
            print(job, st["status"], st["elapsed"], st.get("log_tail", "")[-3:] if st["status"] != "ok" else "")


def compare():
    ref = harness.load(REF)
    out = {}
    for n in CASES:
        if not ((WORK / "A0" / f"{n}.npz").exists() and (WORK / "FIX" / f"{n}.npz").exists()):
            print(n, "skipped (no output: run failed in A0 and/or FIX)")
            continue
        a = harness.load(WORK / "A0" / f"{n}.npz")
        b = harness.load(WORK / "FIX" / f"{n}.npz")
        cab = harness.compare(a, b)
        cra = harness.compare(ref, a)
        crb = harness.compare(ref, b)
        # Omega(a) on the EFT grid vs exact a^2 inside the node range
        ag = a["eft_a_grid"]
        xs = CASES[n]
        lo = float(xs["EFTOmegax1"])
        m = (ag >= max(lo, 1e-2)) & (ag <= 1.0)
        exact = ag[m] ** 2
        row = {
            "A0_vs_FIX_worst_cl": cab["worst_cl"], "A0_vs_FIX_worst_pk": cab["worst_pk"],
            "A0_vs_FIX_sigma8": cab["sigma8"],
            "ref_vs_A0_worst_cl": cra["worst_cl"], "ref_vs_FIX_worst_cl": crb["worst_cl"],
            "ref_vs_A0_sigma8": cra["sigma8"], "ref_vs_FIX_sigma8": crb["sigma8"],
            "Omega_A0_max_rel_err_vs_a2(a>=0.01)": float(max(abs(a["eft_EFTOmegaV"][m] - exact) / exact)),
            "Omega_FIX_max_rel_err_vs_a2(a>=0.01)": float(max(abs(b["eft_EFTOmegaV"][m] - exact) / exact)),
        }
        out[n] = row
        print(n, json.dumps(row, default=str))
    (WORK / "cl_impact.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    {"build": build, "run": run, "compare": compare}[sys.argv[1]](*sys.argv[2:])
