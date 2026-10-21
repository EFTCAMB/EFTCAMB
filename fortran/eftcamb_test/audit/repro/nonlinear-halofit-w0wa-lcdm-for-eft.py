"""Repro (track t8c-features): Halofit / HMcode are fed w0 = -1, wa = 0 for every EFT model with a non-LCDM
background (all models except quintessence with EFTCAMB_effective_w0wa = T).

halofit.f90 EFTCAMB MODs call Params%EFTCAMB%Effective_w_wa, which returns (-1, 0) unless
EFTCAMB_effective_w0wa .and. model%effective_w0wa (09_EFTCAMB_main.f90:700; only 008p4_quintessence sets
effective_w0wa = .true.).  The non-linear correction therefore uses a cosmological-constant expansion history.

Test: pure-EFT designer wCDM with w = -0.8 and every EFT function zero (2_PEFT_wDE_wCDM_1: no modified gravity, DE
perturbations with c_s^2 = 1) vs CAMB's own fluid with w = -0.8, cs2 = 1 (EFTflag = 0).  The two linear P(k) agree to
~1e-5, so the non-linear ratio P_nl/P_lin must agree as well.
Expected: |ratio_EFT / ratio_fluid - 1| ~ 1e-5.  Observed: several per cent at z = 0.
Run: /Users/marco/Software/miniforge3/bin/python nonlinear-halofit-w0wa-lcdm-for-eft.py [TREE]
(TREE defaults to the frozen A0 build; pass a tree with a candidate fix to compare.)
"""

import os
import sys
import tempfile
from contextlib import nullcontext
from pathlib import Path

import numpy as np

os.environ["OMP_NUM_THREADS"] = "2"
TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0")
PARAMS = Path(__file__).resolve().parents[2] / "parameters"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
try:  # hold machine-wide CPU slots of the audit harness when available
    from harness import slot
except Exception:  # noqa: BLE001
    def slot(n):
        return nullcontext()
sys.path.insert(0, TREE)
import camb  # noqa: E402
from camb import model  # noqa: E402

tmp = Path(tempfile.mkdtemp())
COMMON = ("get_tensor_cls = F\nget_scalar_cls = F\nnum_redshiftwindows = 0\nfeedback_level = 0\n"
          "EFTCAMB_write_background = F\n")


def ratio(name, extra, version):
    ini = tmp / f"{name}.ini"
    ini.write_text(f"DEFAULT({PARAMS / name}.ini)\n{COMMON}{extra}")
    cp = camb.read_ini(str(ini))
    cp.set_matter_power(redshifts=[1.0, 0.0], kmax=10.0, nonlinear=True)
    cp.NonLinear = model.NonLinear_pk
    cp.NonLinearModel.set_params(halofit_version=version)
    r = camb.get_results(cp)
    k, z, plin = r.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    _, _, pnl = r.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    return cp, k, z, plin, pnl / plin


print(f"tree: {TREE}")
for version in ("takahashi", "mead2020"):
    with slot(2):
        cpe, k, z, plin_e, r_e = ratio("2_PEFT_wDE_wCDM_1", "", version)
        cpg, k2, _, plin_g, r_g = ratio("1_EFT_GR", "w = -0.8\ncs2_lam = 1\n", version)
    assert np.allclose(k, k2)
    sel = (k > 0.05) & (k < 5)
    print(f"[{version}] EFTflag={cpe.EFTCAMB.EFTflag} ({cpe.EFTCAMB.model_name()}, EFTw0=-0.8) vs "
          f"EFTflag={cpg.EFTCAMB.EFTflag} fluid w={cpg.DarkEnergy.w}")
    for i, zz in enumerate(z):
        dl = np.max(np.abs(plin_e[i, sel] / plin_g[i, sel] - 1))
        dr = r_e[i, sel] / r_g[i, sel] - 1
        j = np.argmax(np.abs(dr))
        print(f"   z={zz:.1f}: linear P(k) max rel diff = {dl:.2e} | P_nl/P_lin EFT vs fluid: max rel diff = "
              f"{dr[j]:+.3e} at k={k[sel][j]:.2f}/Mpc (expected ~1e-5)")
