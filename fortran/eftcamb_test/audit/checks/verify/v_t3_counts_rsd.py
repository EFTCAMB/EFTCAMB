"""Independent check of t3-counts-rsd-gr-branch-missing-de: same tree (A0), EFTflag=0 vs Pure EFT zero functions,
different window (counts z=1.5 sigma 0.05, only the RSD source term on), plus z=0.5 case with all count terms but RSD off.
Also measures RSD-only window with EFT forced on."""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
A0 = AUDIT / "trees" / "A0"
OUT = AUDIT / "work" / "verify" / "t3_counts_rsd"
PE = REPO / "fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa

def worker(tree, ini, out):
    sys.path.insert(0, tree)
    import camb
    cp = camb.read_ini(ini, no_validate=True)
    r = camb.get_results(cp)
    d = r.get_source_cls_dict(raw_cl=True, lmax=600)
    np.savez(out, **{k: np.asarray(v) for k, v in d.items()})

RSDONLY = {"counts_density": "F", "counts_radial": "F", "counts_timedelay": "F", "counts_ISW": "F",
           "counts_velocity": "F", "counts_potential": "F", "counts_redshift": "T",
           "redshift(1)": "1.5", "redshift_sigma(1)": "0.05", "num_redshiftwindows": "1"}

def run(name, ov):
    OUT.mkdir(parents=True, exist_ok=True)
    ini = OUT / f"{name}.ini"
    ini.write_text("\n".join([f"DEFAULT({PE})", f"output_root = {OUT/name}", "EFTCAMB_write_background = F",
                              "feedback_level = 0", "get_tensor_cls = F", "PureEFTmodelOmega = 0"]
                             + [f"{k} = {v}" for k, v in ov.items()]))
    env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "_w", str(A0), str(ini), str(OUT/f"{name}.npz")],
                           env=env, capture_output=True, text=True, cwd=OUT)
    if p.returncode: print(name, "failed", p.stderr[-600:])

def main():
    cases = {"gr_rsd": {"EFTflag": "0", **RSDONLY}, "eft_rsd": RSDONLY,
             "forced_rsd": {**RSDONLY, "EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"},
             "gr_all": {"EFTflag": "0"}, "eft_all": {}}
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda c: run(*c), cases.items()))
    ld = lambda n: dict(np.load(OUT/f"{n}.npz"))
    for a, b in (("eft_rsd", "gr_rsd"), ("forced_rsd", "gr_rsd"), ("eft_all", "gr_all")):
        A, B = ld(a), ld(b)
        for key in [k for k in B if "W1" in k]:
            print(f"{a} vs {b} {key}: max rel {np.max(np.abs(A[key][2:]-B[key][2:]))/np.max(np.abs(B[key][2:])):.3e}")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_w": worker(*sys.argv[2:5])
    else: main()
