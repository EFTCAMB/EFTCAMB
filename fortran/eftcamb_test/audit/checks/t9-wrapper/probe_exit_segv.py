"""Find what makes the process segfault at exit after get_evolution / get_time_evolution calls."""

import sys

from t9lib import WORK, run_snippet

import os
STEM = os.environ.get("PROBE_STEM", "2_PEFT_Omega_power_1")
BASE = f"""
d = WORK / 'dbg' / {STEM!r}
pars = camb.read_ini(str(d / 'wrapper.ini'))
results = camb.get_results(pars)
taus = np.linspace(1.0, 13000.0, 400)
"""
VARIANTS = {
    "none": "",
    "camb_evol": "results.get_time_evolution(1.3, taus, ['delta_cdm'])",
    "eft_evol": "pars.EFTCAMB.get_evolution(results, 1.3, taus)",
    "eft_then_camb": "pars.EFTCAMB.get_evolution(results, 1.3, taus)\nresults.get_time_evolution(1.3, taus, ['delta_cdm'])",
    "camb_then_eft": "results.get_time_evolution(1.3, taus, ['delta_cdm'])\npars.EFTCAMB.get_evolution(results, 1.3, taus)",
    "camb_evol_lowk": "results.get_time_evolution(0.05, taus, ['delta_cdm'])",
    "camb_evol_multi": "results.get_time_evolution([0.01, 0.05, 1.3], taus, ['delta_cdm'])",
}

if __name__ == "__main__":
    tree = sys.argv[1] if len(sys.argv) > 1 else "A0"
    sel = sys.argv[2:] or list(VARIANTS)
    for v in sel:
        r = run_snippet(BASE + VARIANTS[v] + "\nemit(reached_end=True)\n", tree_name=tree, threads=2, timeout=900,
                        cwd=str(WORK / "dbg"))
        print(f"{tree} {v:16s} status={r['status']:16s} reached_end={r['result'].get('reached_end')}")
