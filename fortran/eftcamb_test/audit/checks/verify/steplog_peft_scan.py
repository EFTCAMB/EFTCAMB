"""Adversarial verification (steplog-d4-and-integral): Pure EFT, Omega = steplog(0, v2, aT, delta), small grid.

For each point: run on the fix tree and on A0 (threads=1, <= 4 concurrent). Classifies each point as
  FIX ok & A0 rejected with 'EFTOmegaPPPP is NaN' first  -> falsely rejected by the d4 bug
  both rejected for a physical instability               -> not informative
  both ok                                                -> compare Cls (d4 enters only stability/positivity)
argv[1]: mass stability flag (T/F, default F).
"""

import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "steplog_phys" / "scan"
FIX = harness.AUDIT_ROOT / "work" / "verify" / "steplog_phys" / "fixtree"
INI = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
MASS = sys.argv[1] if len(sys.argv) > 1 else "F"
GRID = list(product([0.02, 0.1, 0.5, 1.0], [0.3, 0.6], [0.3, 0.6]))


def job(args):
    (v2, at, de), tag, tree = args
    name = f"m{MASS}_v{v2}_a{at}_d{de}_{tag}"
    ov = {"EFT_mass_stability": MASS, "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0", "EFTOmega_v2": str(v2),
          "EFTOmega_at": str(at), "EFTOmega_delta": str(de)}
    st = harness.run_model(tree, INI, WORK / f"{name}.npz", overrides=ov, threads=1, timeout=1500)
    log = Path(st["log"]).read_text(errors="replace").splitlines() if st.get("log") else []
    first_nan = next((t.strip() for t in log if "is NaN" in t), "")
    why = [t.strip() for t in log if "Physical instability" in t or "Instability detected" in t][:2]
    return (v2, at, de), tag, {"status": st["status"], "first_nan": first_nan, "why": why}


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    jobs = [(p, t, tr) for p in GRID for t, tr in (("A0", "A0"), ("FIX", FIX))]
    res = {}
    with ThreadPoolExecutor(4) as ex:
        for p, tag, r in ex.map(job, jobs):
            res.setdefault(p, {})[tag] = r
    counts = {"false_reject": 0, "both_unstable": 0, "both_ok": 0, "other": 0}
    rows = []
    for p in GRID:
        a, f = res[p]["A0"], res[p]["FIX"]
        onset = p[1] * math.exp(p[2] ** 2)
        if f["status"] == "ok" and a["status"] != "ok" and a["first_nan"].startswith("EFTOmegaPPPP"):
            cls = "false_reject"
        elif f["status"] != "ok" and a["status"] != "ok":
            cls = "both_unstable"
        elif f["status"] == "ok" and a["status"] == "ok":
            cls = "both_ok"
        else:
            cls = "other"
        counts[cls] += 1
        rows.append({"v2": p[0], "aT": p[1], "delta": p[2], "nan_onset": round(onset, 4), "class": cls,
                     "A0": a, "FIX": f})
        print(f"v2={p[0]:5} aT={p[1]} delta={p[2]} onset={onset:.3f} -> {cls:14s} A0:{a['status']} "
              f"[{a['first_nan'] or ';'.join(a['why'])}] FIX:{f['status']} [{';'.join(f['why'])}]")
    print(counts)
    (WORK.parent / f"scan_mass{MASS}.json").write_text(json.dumps({"counts": counts, "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
