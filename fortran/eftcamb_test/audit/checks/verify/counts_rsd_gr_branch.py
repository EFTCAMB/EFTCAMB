"""Independent check of t3-counts-rsd-gr-branch-missing-de.

Python-API windows (Gaussian counts windows at several z, set programmatically), same-tree comparison
(EFTflag=0 vs Pure EFT with Omega=0) so no B0-vs-A0 differences enter; RSD on/off. Also run on B1 (pre-rebase)
to attribute. If the defect is the missing -2/3 grho_DE term, the RSD-on error should track Omega_DE(z) of the window.
"""

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
OUT = AUDIT / "work" / "verify" / "counts_rsd"
PAR = REPO / "fortran" / "eftcamb_test" / "parameters"
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

ZS = (0.3, 1.0, 2.0)


def worker(tree, eft, rsd, out):
    sys.path.insert(0, tree)
    import camb
    from camb.sources import GaussianSourceWindow

    ini = PAR / ("2_PEFT_Omega_const_1.ini" if eft else "1_EFT_GR.ini")
    cp = camb.read_ini(str(ini), no_validate=True)
    keys = harness.parse_ini(str(ini))
    if eft:
        keys["PureEFTmodelOmega"] = "0"
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    else:
        cp.EFTCAMB.EFTflag = 0
    if Path(tree).name != "B1":
        cp.Reion.max_redshift = 50.0  # same in both arms; B1 cannot reach tau with 50
    cp.WantTensors = False
    cp.InitPower.r = 0.0
    cp.SourceTerms.counts_redshift = bool(rsd)
    cp.SourceWindows = [GaussianSourceWindow(redshift=z, source_type="counts", sigma=0.1, bias=1.0) for z in ZS]
    r = camb.get_results(cp)
    d = r.get_source_cls_dict(raw_cl=True, lmax=600)
    bg = r.get_background_redshift_evolution(np.array(ZS), ["x_e"], format="array")  # noqa: F841
    oml = [float(r.get_Omega("de", z=z)) if hasattr(r, "get_Omega") else np.nan for z in ZS]
    np.savez(out, oml=np.array(oml), eftflag=int(cp.EFTCAMB.EFTflag),
             **{k: np.asarray(v) for k, v in d.items()})


def run(tree, eft, rsd):
    name = f"{Path(tree).name}_{'eft' if eft else 'gr'}_{'rsd' if rsd else 'nor'}"
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "_worker", str(tree), str(int(eft)), str(int(rsd)),
                            str(OUT / f"{name}.npz")], env=env, capture_output=True, text=True, cwd=OUT)
    if p.returncode:
        print(name, "failed", p.stderr[-1500:])
    return name


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trees = [Path(a) if "/" in a else AUDIT / "trees" / a for a in (sys.argv[1:] or ["A0", "B1"])]
    cases = [(t, e, r) for t in trees for e in (0, 1) for r in (0, 1)]
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda c: run(*c), cases))
    for t in trees:
        for rsd in (0, 1):
            tag = "rsd" if rsd else "nor"
            fa, fb = OUT / f"{t.name}_eft_{tag}.npz", OUT / f"{t.name}_gr_{tag}.npz"
            if not (fa.exists() and fb.exists()):
                continue
            a, b = np.load(fa), np.load(fb)
            row = {}
            for i in range(1, len(ZS) + 1):
                k = f"W{i}xW{i}"
                row[f"z={ZS[i-1]}"] = float(np.max(np.abs(a[k][2:] - b[k][2:]) / np.abs(b[k][2:])))
            print(f"{t.name:4s} {tag}: eftflag={int(a['eftflag'])} OmegaDE(z)={np.round(b['oml'], 3).tolist()} "
                  + "  ".join(f"{k}:{v:.2e}" for k, v in row.items()))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    else:
        main()
