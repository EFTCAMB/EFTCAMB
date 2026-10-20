"""t3-gr (a): localize residual rounding-level differences (variant tree vs B0) by switching ini features off."""
import sys
from concurrent.futures import ThreadPoolExecutor
from common import harness, run, fmt, WORK

name = sys.argv[1]
tree = WORK / name if name not in harness.TREES else name
sub = "a4"
cases = {
    "nonu": {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0"},
    "nonu_notens": {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0",
                    "get_tensor_cls": "F"},
    "nolens": {"do_lensing": "F"},
    "notrans": {"get_transfer": "F"},
    "noreion": {"reionization": "F"},
}
if len(sys.argv) > 2:
    cases = {k: cases[k] for k in sys.argv[2].split(",")}
jobs = []
with ThreadPoolExecutor(4) as ex:
    for c, ov in cases.items():
        ps = {"Reion.max_redshift": 50.0}
        jobs.append((c, ex.submit(run, "B0", f"B0_{c}", overrides=ov, sub=sub),
                     ex.submit(run, tree, f"{name}_{c}", overrides=ov, python_set=ps, sub=sub)))
    for c, fb, fv in jobs:
        r = harness.compare(fv.result(), fb.result())
        print(f"{c:12s}", fmt(r), "| unl:", {s: f"{x:.0e}" for s, x in r['cl'].get('unlensed_scalar', {}).items()},
              "PP:", f"{r['cl'].get('lens_potential', {}).get('PP', 0):.0e}")
