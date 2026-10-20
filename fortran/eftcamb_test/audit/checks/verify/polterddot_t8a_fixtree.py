"""Verifier (t8a polterddot-nonflat-gr): attribution through the finding's proposed fix, alternative 2.

Builds a copy of A0 in which polterddot is CAMB's original expression whenever
`CP%EFTCAMB%EFTflag == 0 .or. .not. EV%EFTCAMBactive` (EFTCAMB's replacement kept for active EFT), then runs the
1_EFT_GR.ini test ini (ini code path, massive nu, tau=0.09, lensing on; different from the repro's set_cosmology run)
for several omk on A0, B0 and FIX, all with Reion.max_redshift=50 (removes the separately known reionization default
change). Also checks two flat EFT models are unchanged by the fix (round-off only).
Usage: python polterddot_t8a_fixtree.py build | run
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "polterddot_t8a"
FIX = WORK / "fixtree"

OLD = (
    "            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*octgdot - dopacity*(pig - polter) - &\n"
    "                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot\n"
    "            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*Edot(3)\n"
    "            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot\n"
)
NEW = (
    "            if (CP%EFTCAMB%EFTflag == 0 .or. .not. EV%EFTCAMBactive) then\n"
    "                polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &\n"
    "                    k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &\n"
    "                    50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &\n"
    "                    *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &\n"
    "                    polter - 3.0d0/5.0d0*E(2))*dopacity + (-1.0d0/10.0d0*pigdot &\n"
    "                    + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity\n"
    "            else\n"
    "            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*octgdot - dopacity*(pig - polter) - &\n"
    "                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot\n"
    "            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*Edot(3)\n"
    "            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot\n"
    "            end if\n"
)


def build():
    import shutil

    if FIX.exists():
        shutil.rmtree(FIX)
    shutil.copytree(harness.tree_path("A0"), FIX, symlinks=True)
    f = FIX / "fortran" / "equations.f90"
    s = f.read_text()
    assert s.count(OLD) == 1
    f.write_text(s.replace(OLD, NEW))
    print(json.dumps(harness.build_tree(FIX, FIX), indent=1))


GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OV = {"get_tensor_cls": "F", "l_max_scalar": "2000", "get_transfer": "F", "num_redshiftwindows": "0"}
OMKS = ("0", "-0.05", "0.05", "-0.1", "0.1")
EFT_FLAT = ("2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1")


def job(a):
    tree, ini, omk, name = a
    ov = dict(OV)
    if omk is not None:
        ov["omk"] = omk
    return harness.run_model(tree, ini, WORK / "runs" / f"{name}.npz", overrides=ov,
                             python_set={"Reion.max_redshift": 50.0}, threads=1)


def stats(x, b, key):
    a, b = x[key], b[key]
    n = min(len(a), len(b), 2001)
    tt = np.abs(a[2:n, 0] / b[2:n, 0] - 1)
    te = np.abs(a[2:n, 3] - b[2:n, 3]) / np.sqrt(b[2:n, 0] * b[2:n, 1])
    ee = np.abs(a[2:n, 1] / b[2:n, 1] - 1)
    return dict(TT_l2=float(tt[0]), TT_max=float(tt.max()), TT_argmax=int(2 + np.argmax(tt)),
                TT_max_l_ge100=float(tt[98:].max()), TE_max=float(te.max()), EE_max=float(ee.max()))


def run():
    trees = {"A0": "A0", "B0": "B0", "FIX": str(FIX)}
    jobs = [(t, GR, o, f"{lab}_GR_omk{o}") for o in OMKS for lab, t in trees.items()]
    jobs += [(t, harness.PARAMS_DIR / f"{m}.ini", None, f"{lab}_{m}") for m in EFT_FLAT
             for lab, t in trees.items() if lab != "B0"]
    with ThreadPoolExecutor(4) as ex:
        sts = list(ex.map(job, jobs))
    for st in sts:
        if st["status"] != "ok":
            print("FAILED", st)
    L = lambda n: harness.load(WORK / "runs" / f"{n}.npz")  # noqa: E731
    out = {}
    for o in OMKS:
        b0 = L(f"B0_GR_omk{o}")
        for lab in ("A0", "FIX"):
            x = L(f"{lab}_GR_omk{o}")
            for key in ("cl_unlensed_scalar", "cl_lensed_scalar"):
                s = stats(x, b0, key)
                out[f"omk={o} {lab} vs B0 {key}"] = s
                print(f"omk={o:>6s} {lab:>3s} vs B0 {key:>18s}: TT(l=2) {s['TT_l2']:.3e} max TT {s['TT_max']:.3e}"
                      f" (l={s['TT_argmax']}) TT l>=100 {s['TT_max_l_ge100']:.3e} TE/sqrt {s['TE_max']:.3e}"
                      f" EE {s['EE_max']:.3e}")
    for m in EFT_FLAT:
        c = harness.compare(L(f"FIX_{m}"), L(f"A0_{m}"), lmax=2000)
        out[f"{m} FIX vs A0"] = {"worst_cl": c["worst_cl"], "bit_identical": c["bit_identical"]}
        print(f"{m}: FIX vs A0 worst_cl {c['worst_cl']} bit_identical {c['bit_identical']}")
    (WORK / "fixtree_summary.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    {"build": build, "run": run}[sys.argv[1]]()
