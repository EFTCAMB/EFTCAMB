"""t3-gr (a): EFTflag=0 code under audit vs plain CAMB v2 over a set of GR configurations.

python a5_variants.py TREE_A TREE_B [case,case,...] [--maxz50]
TREE_A/TREE_B: harness tree names or variant names in the t3-gr work dir.
--maxz50: set Reion.max_redshift=50 on both sides (removes the known reionization-bracket change).
Writes a json summary to $WORK/a5_<A>_vs_<B>[_maxz50].json
"""

from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor

from common import GR_INI, INIS, WORK, dump, fmt, harness, run, summary

NONU = {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0"}
CASES: dict[str, dict[str, str]] = {
    "base": {},
    "nonu": NONU,
    "nu3": {"omnuh2": "0.0064", "massive_neutrinos": "3", "massless_neutrinos": "0.046", "nu_mass_eigenstates": "1"},
    "omk_p05": {"omk": "0.05", "num_redshiftwindows": "0"},
    "omk_m05": {"omk": "-0.05", "num_redshiftwindows": "0"},
    "omk_m01_nonu": {"omk": "-0.01", "num_redshiftwindows": "0", **NONU},
    "tensors_only": {"get_scalar_cls": "F", "get_transfer": "F", "do_lensing": "F"},
    "nolens": {"do_lensing": "F"},
    "lens_nl": {"do_nonlinear": "2", "halofit_version": "9"},
    "pk_nl": {"do_nonlinear": "3", "halofit_version": "4"},
    "multi_z": {"transfer_num_redshifts": "4", "transfer_redshift(1)": "3", "transfer_redshift(2)": "1",
                "transfer_redshift(3)": "0.5", "transfer_redshift(4)": "0",
                "transfer_filename(2)": "t2.dat", "transfer_filename(3)": "t3.dat", "transfer_filename(4)": "t4.dat",
                "transfer_matterpower(2)": "p2.dat", "transfer_matterpower(3)": "p3.dat",
                "transfer_matterpower(4)": "p4.dat"},
    "kmax_small": {"transfer_kmax": "0.01", "do_lensing": "F", "get_tensor_cls": "F"},
    "kmax_small_pk": {"transfer_kmax": "0.01", "do_lensing": "F", "get_tensor_cls": "F", "get_scalar_cls": "F",
                      "num_redshiftwindows": "0"},
    "w09_fluid": {"w": "-0.9", "dark_energy_model": "fluid"},
    "w11_fluid": {"w": "-1.1", "dark_energy_model": "ppf"},
    "w09_wa_ppf": {"w": "-0.9", "wa": "-0.3", "dark_energy_model": "ppf"},
    "acc2": {"accuracy_boost": "2", "l_accuracy_boost": "2"},
    "zre_input": {"re_use_optical_depth": "F", "re_redshift": "9"},
    "theta_ini": {"__ini__": str(INIS / "t3gr_GR_theta.ini")},
}


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    maxz = "--maxz50" in sys.argv
    ta, tb = args[0], args[1]
    cases = list(CASES) if len(args) < 3 else args[2].split(",")
    tree = {t: (t if t in harness.TREES else WORK / t) for t in (ta, tb)}
    ps = {"Reion.max_redshift": 50.0} if maxz else None
    sub = "a5" + ("_maxz50" if maxz else "")
    res = {}
    with ThreadPoolExecutor(int(__import__("os").environ.get("T3_JOBS", "4"))) as ex:
        futs = {}
        for c in cases:
            ov = {k: v for k, v in CASES[c].items() if k != "__ini__"}
            ini = CASES[c].get("__ini__", GR_INI)
            for t in (ta, tb):
                futs[(c, t)] = ex.submit(run, tree[t], f"{t}_{c}", ini=ini, overrides=ov, python_set=ps, sub=sub)
        for c in cases:
            pa, pb = futs[(c, ta)].result(), futs[(c, tb)].result()
            try:
                r = harness.compare(pa, pb)
            except Exception as e:  # noqa: BLE001
                print(f"{c:14s} compare failed: {e!r}")
                continue
            res[c] = summary(r)
            print(f"{c:14s}", fmt(r))
    dump(res, WORK / f"a5_{ta}_vs_{tb}{'_maxz50' if maxz else ''}.json")


if __name__ == "__main__":
    main()
