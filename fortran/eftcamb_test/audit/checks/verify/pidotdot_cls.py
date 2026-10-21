"""CMB-level verification of output-pidotdot-zero through the audit harness (lensed + unlensed Cls, P(k)).

A0 (code under audit) vs FX (A0 + proposed fix: Einstein factors and dz refreshed after the pi equation).
Affected: Gamma4 + 2 Gamma5 != 0 with pi != 0 (Horava, the verify PEFT Omega+Gamma5 model), and Omega'/Gamma2 != 0
via EFTeomM -> dz (f(R)). Controls: GR, Horndeski-type models.

Run: /Users/marco/Software/miniforge3/bin/python pidotdot_cls.py [chunk-index]
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

AUDIT = harness.AUDIT_ROOT
FX = AUDIT / "work/t8a-equations/FX"
OUT = AUDIT / "work/verify/pidotdot/cls"
P = harness.PARAMS_DIR
I = HERE.parents[1] / "inis"  # noqa: E741

CHUNKS = [
    [P / f"5_Horava_eta_{i}.ini" for i in (1, 2, 3, 4)] + [I / "verify_pidd_omega_g5m.ini", I / "verify_pidd_omega_only.ini"],
    [P / f"5_Horava_SolSyst_{i}.ini" for i in (1, 2, 3)] + [P / "4_Des_F_R_LCDM_1.ini", P / "1_EFT_GR.ini"],
    [P / "3_RPH_KT_power_1.ini", P / "5_hdsk_TG_2.ini", P / "2_PEFT_gamma4_const_1.ini", P / "2_PEFT_gamma5_const_1.ini",
     P / "2_PEFT_Omega_lin_1.ini", I / "008p5-beyond-horndeski.ini"],
]


def cl_report(a, b):
    """fractional differences of unlensed scalar TT (b = reference)."""
    n = min(len(a["cl_unlensed_scalar"]), len(b["cl_unlensed_scalar"]), 2501)
    ca, cb = a["cl_unlensed_scalar"][:n], b["cl_unlensed_scalar"][:n]
    ell = np.arange(n)
    m = ell >= 2
    dtt = (ca[m, 0] - cb[m, 0]) / cb[m, 0]
    i = int(np.argmax(np.abs(dtt)))
    pick = {L: float(dtt[L - 2]) for L in (2, 3, 4, 10, 30, 100, 1000) if L < n}
    return {"max_dTT": float(np.max(np.abs(dtt))), "ell_at_max": int(ell[m][i]), "dTT_at": pick}


def main():
    idx = [int(x) for x in sys.argv[1:]] or range(len(CHUNKS))
    OUT.mkdir(parents=True, exist_ok=True)
    report = {}
    for c in idx:
        inis = CHUNKS[c]
        need_a0 = [p for p in inis if not (AUDIT / "runs/A0" / f"{p.stem}.npz").exists()]
        todo_fx = [p for p in inis if not (OUT / f"{p.stem}_FX.npz").exists()]
        need_a0 = [p for p in need_a0 if not (OUT / f"{p.stem}_A0.npz").exists()]
        st_fx = harness.run_batch(FX, todo_fx, OUT, jobs=4, suffix="_FX") if todo_fx else []
        st_a0 = harness.run_batch("A0", need_a0, OUT, jobs=4, suffix="_A0") if need_a0 else []
        status = {Path(s["model_ini"]).stem + "_FX": s["status"] for s in st_fx}
        status.update({Path(s["model_ini"]).stem + "_A0": s["status"] for s in st_a0})
        for p in inis:
            ref = AUDIT / "runs/A0" / f"{p.stem}.npz"
            if not ref.exists():
                ref = OUT / f"{p.stem}_A0.npz"
            fx = OUT / f"{p.stem}_FX.npz"
            if not (ref.exists() and fx.exists()):
                report[p.stem] = {"error": status}
                print(f"{p.stem:32s} ERROR {status}")
                continue
            a, f = harness.load(ref), harness.load(fx)
            cmp = harness.compare(f, a)
            r = cl_report(f, a)
            r.update({"worst_cl": cmp["worst_cl"], "worst_pk": cmp["worst_pk"], "bit_identical": cmp["bit_identical"],
                      "lensed_TT": cmp["cl"].get("total", {}).get("TT"), "EE": cmp["cl"].get("total", {}).get("EE"),
                      "ref": str(ref)})
            report[p.stem] = r
            print(f"{p.stem:32s} FX vs A0: unlensed max|dTT/TT|={r['max_dTT']:.2e} (ell={r['ell_at_max']}) "
                  f"dTT(2,4,30,1000)={r['dTT_at'].get(2, 0):+.1e},{r['dTT_at'].get(4, 0):+.1e},"
                  f"{r['dTT_at'].get(30, 0):+.1e},{r['dTT_at'].get(1000, 0):+.1e} "
                  f"lensedEE={r['EE'] if r['EE'] is None else format(r['EE'], '.1e')} "
                  f"worst_pk={r['worst_pk']} bit_identical={r['bit_identical']}")
        (OUT / f"report_chunk{'_'.join(map(str, idx))}.json").write_text(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()
