"""Test the proposed fix for 'spline-derivatives-outside-nodes' in a private copy of A0 (verifier's own patch).

Patch: in 04p016 (Spline1D{First..Fourth}Derivative) and 04p019 (SplineQ1D{First..Fourth}Derivative) the result is
initialised to 0._dl instead of self%null_value; the value functions keep null_value.  Then rerun the partial-coverage
cases of verify_spline_outside_nodes.py on the fixed tree and compare with EFTwDE = 0.

usage: python verify_spline_outside_fix.py [build|run|analyse]   (default: all)
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness  # noqa: E402
import verify_spline_outside_nodes as V  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "spline_outside_fix"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "spline_outside_fixruns"
SRC = "fortran/eftcamb/04f_parametrizations_1D/"
PATCH = {
    "04p016_spline_parametrizations_1D.f90": [f"Spline1D{d}Derivative = self%null_value"
                                              for d in ("First", "Second", "Third", "Fourth")],
    "04p019_spline5_parametrizations_1D.f90": [f"SplineQ1D{d}Derivative = self%null_value"
                                               for d in ("First", "Second", "Third", "Fourth")],
}
CASES = ["ref", "s9_t4a", "s10_t4a", "s10_lategap05", "s10_earlygap", "s9_endpoint1"]


def build():
    if not DEST.exists():
        shutil.copytree(harness.TREES["A0"], DEST, symlinks=True)
        for fname, olds in PATCH.items():
            f = DEST / SRC / fname
            s = f.read_text()
            for old in olds:
                assert s.count(old) == 1, (fname, old, s.count(old))
                s = s.replace(old, old.replace("self%null_value", "0._dl"))
            f.write_text(s)
    print(harness.build_tree("A0", DEST))


def run():
    from concurrent.futures import ThreadPoolExecutor

    OUT.mkdir(parents=True, exist_ok=True)

    def one(name):
        return name, harness.run_model(DEST, V.INI, OUT / f"{name}.npz", overrides=V.CASES[name], threads=1)

    todo = [c for c in CASES if not (OUT / f"{c}.npz").exists()]
    with ThreadPoolExecutor(4) as ex:
        for name, st in ex.map(one, todo):
            print(f"{name:16s} {st['status']}  {st['elapsed']}s")


def analyse():
    refA0 = harness.load(V.OUT / "ref.npz")
    c = harness.compare(harness.load(OUT / "ref.npz"), refA0)
    print(f"FIX ref (EFTwDE=0) vs A0 ref: bit_identical={c['bit_identical']} worst_cl={c['worst_cl']:.2e}")
    print(f"{'case':16s} {'A0 worstCl':>11s} {'FIX worstCl':>11s} {'A0 P(k)':>10s} {'FIX P(k)':>10s} "
          f"{'FIX Ldot':>10s} FIX bit_identical")
    for name in CASES[1:]:
        a0 = harness.compare(harness.load(V.OUT / f"{name}.npz"), refA0)
        fx = harness.compare(harness.load(OUT / f"{name}.npz"), refA0)
        print(f"{name:16s} {a0['worst_cl']:11.3e} {fx['worst_cl']:11.3e} {a0['worst_pk']:10.3e} {fx['worst_pk']:10.3e} "
              f"{fx['eft']['EFTLambdadot']:10.3e} {fx['bit_identical']}")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("build", "all"):
        build()
    if what in ("run", "all"):
        run()
    if what in ("analyse", "all"):
        analyse()
