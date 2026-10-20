"""Independent check of eft-turnon-jump-cmb-time-sampling: sweep EFTCAMB_turn_on_time finely around 0.01.
Physical dependence on a_on is smooth; a first-order quadrature error of a jump on a coarse grid gives a
non-smooth (sawtooth) residual whose amplitude shrinks with TimeStepBoost."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "turnon_sweep"
OUT.mkdir(parents=True, exist_ok=True)
MODEL = sys.argv[1] if len(sys.argv) > 1 else "2_PEFT_Omega_const_1"
AON = np.linspace(0.0096, 0.0104, 9)
BOOSTS = ["1", "4"]
ELLS = [2, 5, 10, 30]


def job(args):
    aon, tsb = args
    out = OUT / f"{MODEL}_a{aon:.6f}_ts{tsb}.npz"
    over = {"EFTCAMB_turn_on_time": f"{aon:.8f}", "TimeStepBoost": tsb, "num_redshiftwindows": "0",
            "l_max_scalar": "200", "get_tensor_cls": "F", "do_lensing": "F", "get_transfer": "F"}
    if not out.exists():
        st = harness.run_model("A0", harness.PARAMS_DIR / f"{MODEL}.ini", out, overrides=over, threads=1, timeout=1200)
        if st["status"] != "ok":
            raise SystemExit(f"{aon} {tsb}: {st['status']} {st.get('log_tail')}")
    return (aon, tsb), harness.load(out)["cl_unlensed_scalar"][:, 0]


with ThreadPoolExecutor(4) as ex:
    cl = dict(ex.map(job, [(a, b) for b in BOOSTS for a in AON]))
for tsb in BOOSTS:
    print(f"{MODEL} TimeStepBoost={tsb}")
    for l in ELLS:
        y = np.array([cl[(a, tsb)][l] for a in AON])
        # residual from a quadratic fit = non-smooth part
        p = np.polyfit(AON - 0.01, y, 2)
        res = y / np.polyval(p, AON - 0.01) - 1
        d2 = np.diff(y, 2) / y[1:-1]
        print(f"  ell={l:3d} TT/TT(a=0.01)-1:", " ".join(f"{v/ y[4]-1:+.3e}" for v in y))
        print(f"          quad-fit resid max {np.abs(res).max():.2e}; max|2nd diff| {np.abs(d2).max():.2e}")
for l in ELLS:
    r = [cl[(a, '1')][l] / cl[(a, '4')][l] - 1 for a in AON]
    print(f"ell={l:3d} ts1/ts4-1 vs a_on:", " ".join(f"{v:+.2e}" for v in r))
