"""Uninitialised timestep-cache components: run models on a build with -finit-real=nan -finit-derived (every
uninitialised local real, including derived-type components, starts as NaN) and report the stability verdict and the
NaN-check message (feedback 3)."""
import sys, time
sys.path.insert(0, ".")
import t6lib, harness
NAN = str(t6lib.WORK / "nanini")
models = sys.argv[1:] or ["2_PEFT_gamma1_const_1", "5_quint_exponential"]
for m in models:
    ini = harness.PARAMS_DIR / f"{m}.ini"
    for tree in ("A0", NAN):
        t0 = time.time()
        r = t6lib.probe(tree, ini, overrides={"feedback_level": "3", "EFTCAMB_skip_stability": "F"}, timeout=240)
        msgs = sorted({l.strip() for l in r["log"] if "Nan" in l or "nan" in l.lower() or "instab" in l.lower()
                       or "theory" in l})
        print(f"{m:24s} {('A0' if tree == 'A0' else 'nanini'):7s} status={r['status']} {time.time()-t0:6.1f}s "
              f"verdict={'stable' if r.get('bg_error') is None else 'UNSTABLE'} msgs={msgs[:6]}")
