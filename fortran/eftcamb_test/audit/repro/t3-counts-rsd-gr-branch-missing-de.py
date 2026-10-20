"""Repro (t3-gr): number-counts redshift-distortion source is wrong for EFT models while EFT perturbations are off.

For EFTflag/=0, derivs sets grhov_t = 0 and keeps grho, gpres without the dark-energy part (the background comes from
the EFT model). When EV%EFTCAMBactive is false (before the RGR/turn-on time, or for a model at its GR point, whose
RGR time is 1.1 so EFT never switches on) output_window_sources uses the GR branch
    counts_redshift_source = ((4*adotoa**2 + gpres + grho/3)/k*wing2 + ...)*sigma + ...
which then misses gpres_DE + grho_DE/3 = -(2/3)*grhov_t for a cosmological constant.

Test: Pure EFT with all EFT functions zero (= LCDM) vs plain CAMB v2 (B0) LCDM, base_params.ini number-count window
(z=0.5, counts_redshift=T) and lensing window. Control rows: same with counts_redshift=F, and with the EFT path forced
on (EFTCAMB_GR_threshold=-1, stability check skipped since all EFT functions vanish), where the EFT branch of the
source is used. Optional argument: a tree with the proposed fix (GR-branch RSD term from the Friedmann equations,
4 adotoa^2 + gpres + grho/3 -> 4 adotoa^2 - 2 Hdot when EFTflag/=0).
Expected: all source-window Cls agree with B0 to <~1e-4. Observed: printed.
"""

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
TREES = {"A0": AUDIT / "trees" / "A0", "B0": AUDIT / "trees" / "B0"}
PAR = REPO / "fortran" / "eftcamb_test" / "parameters"
OUT = AUDIT / "work" / "t3-gr" / "repro_counts"
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402  (machine-wide CPU slots)


def worker(tree, ini, out):
    sys.path.insert(0, tree)
    import camb

    cp = camb.read_ini(ini, no_validate=True)
    if hasattr(cp, "Reion"):
        cp.Reion.max_redshift = 50.0  # remove the (separate) reionization-bracket difference
    r = camb.get_results(cp)
    d = r.get_source_cls_dict(raw_cl=True, lmax=800)
    np.savez(out, **{k: np.asarray(v) for k, v in d.items()})


def run(name, tree, base, ov):
    OUT.mkdir(parents=True, exist_ok=True)
    ini = OUT / f"{name}.ini"
    ini.write_text("\n".join([f"DEFAULT({base})", f"output_root = {OUT / name}", "EFTCAMB_write_background = F",
                              "feedback_level = 0", "get_tensor_cls = F"] + [f"{k} = {v}" for k, v in ov.items()]))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "_worker", str(TREES[tree]), str(ini), str(OUT / f"{name}.npz")],
                           env=env, capture_output=True, text=True, cwd=OUT)
    if p.returncode:
        print(name, "failed", p.stderr[-800:])


def main():
    pe = PAR / "2_PEFT_Omega_const_1.ini"
    gr = PAR / "1_EFT_GR.ini"
    cases = [("B0", "B0", gr, {}), ("B0_noRSD", "B0", gr, {"counts_redshift": "F"}),
             ("PEFTzero", "A0", pe, {"PureEFTmodelOmega": "0"}),
             ("PEFTzero_noRSD", "A0", pe, {"PureEFTmodelOmega": "0", "counts_redshift": "F"}),
             ("PEFTzero_forcedEFT", "A0", pe, {"PureEFTmodelOmega": "0", "EFTCAMB_GR_threshold": "-1",
                                                "EFTCAMB_skip_stability": "T"})]
    if len(sys.argv) > 1:  # optional tree with the proposed fix
        TREES["FIX"] = Path(sys.argv[1])
        cases.append(("PEFTzero_fixed", "FIX", pe, {"PureEFTmodelOmega": "0"}))
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda c: run(c[0], c[1], c[2], c[3]), cases))

    def ld(n):
        with np.load(OUT / f"{n}.npz") as f:
            return {k: f[k] for k in f.files}

    for n, ref in (("PEFTzero", "B0"), ("PEFTzero_noRSD", "B0_noRSD"), ("PEFTzero_forcedEFT", "B0"),
                   ("PEFTzero_fixed", "B0")):
        if not (OUT / f"{n}.npz").exists():
            continue
        a, b = ld(n), ld(ref)
        row = {k: float(np.max(np.abs(a[k][2:] - b[k][2:])) / np.max(np.abs(b[k][2:]))) for k in ("W1xW1", "W1xW2",
                                                                                                   "W2xW2", "PxW1")}
        print(f"{n:20s} vs {ref:9s}: " + "  ".join(f"{k}={v:.2e}" for k, v in row.items()))
    print("expected: all <~1e-4 (W1 = number counts z=0.5 with RSD, W2 = lensing window)")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(*sys.argv[2:5])
    else:
        main()
