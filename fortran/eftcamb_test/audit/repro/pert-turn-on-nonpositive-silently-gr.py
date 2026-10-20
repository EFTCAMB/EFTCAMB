"""Repro (t9-wrapper): EFTCAMB_turn_on_time (EFTCAMB_pert_turn_on) = 0, negative, NaN or > 1 is accepted silently and
the EFT perturbations are never switched on: the result is the EFT background with GR perturbations, no error.

Mechanism: unless EFTCAMB_skip_RGR = T, CAMBdata_SetParams replaces the turn-on by the return-to-GR time
(results.f90:599-604, this%CP%EFTCAMB%EFTCAMB_pert_turn_on = RGR_time). EFTCAMBReturnToGR (09_EFTCAMB_RGR.f90:143-171)
scans a log10 grid from log10(a_on) to 0: for a_on = 0 the grid is -inf/NaN, for a_on < 0 or NaN it is NaN, so no
EFT function ever exceeds the GR threshold and RGR_time = 1.1 (after today) becomes the turn-on time; for a_on > 1
the first grid point already exceeds the threshold and RGR_time = a_on > 1.  The perturbation switch then never
fires (EFTCAMBactive stays .false.).  The only validation is the "EFTCAMB starts before thermo tauminn" check in
initial() (equations.f90:2186-2192), which catches tiny positive values (1e-30 below) but not these.
With EFTCAMB_skip_RGR = T the same 0 and -0.1 do raise the tauminn error (last rows): the RGR step hides them.
The python docstring of the sister flag EFTCAMB_back_turn_on documents 0 as "at all times", so 0 is a natural input.

Model: Pure EFT, Omega = 0.01 constant, LCDM background (fresh CAMBparams, lmax 1000).
Expected: an error for a_on <= 0, NaN or > 1 (or 0 meaning "from the start"), never a silent GR-perturbation result.
Observed: printed below; max|C_l^TT/C_l^TT(GR) - 1| ~ 4e-4 (the background-only difference) instead of ~3e-2.

Run: /Users/marco/Software/miniforge3/bin/python pert-turn-on-nonpositive-silently-gr.py [TREE]
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0")
sys.path.insert(0, TREE)
import numpy as np  # noqa: E402

import camb  # noqa: E402


def tt(eft):
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    p.set_for_lmax(1000, lens_potential_accuracy=0)
    p.EFTCAMB.initialize_parameters(p, eft, print_header=False)
    r = camb.get_results(p)
    return r.get_cmb_power_spectra(raw_cl=True, CMB_unit="muK")["total"][2:1000, 0]


base = {"EFTflag": 1, "PureEFTmodel": 1, "PureEFTmodelOmega": 1, "EFTOmega0": 0.01, "EFTwDE": 0}
ref = tt(base)  # default turn-on a = 0.01
gr = tt({"EFTflag": 0})
print(f"camb from {camb.__file__}")
print(f"reference: default turn-on a=0.01; GR vs reference max|dTT/TT| = {np.max(np.abs(gr / ref - 1)):.3e}")
print(f"{'EFTCAMB_turn_on_time':>22s}  {'outcome':40s} {'max|TT/TT_ref-1|':>17s} {'max|TT/TT_GR-1|':>16s}")
for v, skip in ((1e-3, False), (1e-30, False), (0.0, False), (-0.1, False), (float("nan"), False), (2.0, False),
                (0.0, True), (-0.1, True)):
    try:
        c = tt(dict(base, EFTCAMB_turn_on_time=v, EFTCAMB_skip_RGR=skip))
    except camb.CAMBError as e:
        print(f"{v:>16g}{' noRGR' if skip else '      '}  {'CAMBError: ' + str(e).splitlines()[-1][:28]:40s}")
        continue
    print(f"{v:>16g}{' noRGR' if skip else '      '}  {'no error, result returned':40s} {np.max(np.abs(c / ref - 1)):17.3e} "
          f"{np.max(np.abs(c / gr - 1)):16.3e}")
print("expected for 0, -0.1, nan, 2: an error (as for 1e-30), observed: silent result equal to GR perturbations")
