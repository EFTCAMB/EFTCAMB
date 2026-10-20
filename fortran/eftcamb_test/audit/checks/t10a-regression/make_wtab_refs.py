"""Create EFTflag=0 GR references with CAMB's own fluid dark energy (cs2=1) and w(a) tabulated from the A0 EFT
background (eft_gpresq/eft_grhoq) of a model, for models whose effective dark-energy density stays positive with
w >= -1 (so the fluid is well defined).  Used as matched-background references in double_ratio.py.

Usage: python make_wtab_refs.py [--eff] MODEL [MODEL ...]

--eff: use the GR-equivalent dark energy of the model's expansion history, rho_eff = 3 H^2 - rho_m (any Omega),
       with w_eff from continuity, w = -1 - dlog(rho_eff a^-2)/dlog(a)/3; the reference uses the PPF dark energy
       model when w_eff crosses -1 (fluid otherwise).  Files are named t10a_GR_weff_<model>.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

INIS = Path(__file__).resolve().parents[2] / "inis"
RUNS_A0 = Path(harness.AUDIT_ROOT) / "runs" / "A0"


def main_eff(models) -> None:
    from scipy.interpolate import CubicSpline

    for m in models:
        d = harness.load(RUNS_A0 / f"{m}.npz")
        a = d["eft_a_grid"]
        r = 3 * d["eft_adotoa"] ** 2 - d["eft_grhom_t"]
        frac = r / (3 * d["eft_adotoa"] ** 2)
        if np.any(r <= 0):
            print(f"{m}: rho_eff <= 0 somewhere (min fraction {frac.min():.2e}); no reference possible")
            continue
        la = np.log(a)
        w = -1 - CubicSpline(la, np.log(r / a**2))(la, 1) / 3
        # where the GR-equivalent dark energy is negligible the subtraction is noisy: freeze w there
        bad = frac < 1e-5
        if np.any(bad):
            i0 = np.max(np.nonzero(bad)[0]) + 1
            w[:i0] = w[i0]
        w = np.clip(w, -3.0, 1.0)
        model = "ppf" if np.any(w < -1) else "fluid"
        dat = INIS / f"t10a_weff_{m}.dat"
        hdr = f"# a w_eff  (GR-equivalent dark energy of the A0 expansion history of {m}: 3H^2 - rho_m)\n"
        dat.write_text(hdr + "".join(f"{x:.16e} {y:.16e}\n" for x, y in zip(a, w)))
        base = "hdsk_base_params.ini" if "hdsk_base" in (harness.PARAMS_DIR / f"{m}.ini").read_text() else "base_params.ini"
        ini = harness.parse_ini(harness.PARAMS_DIR / f"{m}.ini")
        extra = f"omnuh2 = {ini.get('omnuh2')}\n" if m.startswith(("5_Kmouflage", "5_Kmimic")) else ""
        txt = (f"# t10a: EFTflag=0 reference, {model} dark energy with w_eff(a) of the A0 expansion history of {m}\n"
               f"EFTflag = 0\ndark_energy_model = {model}\nuse_tabulated_w = T\nwafile = {dat}\ncs2_lam = 1\n{extra}"
               f"output_root = t10a_GR_weff_{m}\nDEFAULT({harness.PARAMS_DIR / base})\n")
        (INIS / f"t10a_GR_weff_{m}.ini").write_text(txt)
        print(f"{m}: {model}, w_eff in [{w.min():.3f}, {w.max():.3f}], min DE fraction {frac.min():.1e}")


def main() -> None:
    if sys.argv[1] == "--eff":
        main_eff(sys.argv[2:])
        return
    for m in sys.argv[1:]:
        d = harness.load(RUNS_A0 / f"{m}.npz")
        a, rq, pq = d["eft_a_grid"], d["eft_grhoq"], d["eft_gpresq"]
        if np.max(np.abs(d["eft_EFTOmegaV"])) > 0:
            print(f"{m}: Omega != 0, grhoq is not the GR-equivalent dark energy; no fluid reference")
            continue
        if np.any(rq <= 0):
            print(f"{m}: grhoq <= 0 somewhere, fluid reference not possible")
            continue
        w = pq / rq
        if np.any(w < -1 - 1e-12):
            print(f"{m}: w < -1 somewhere (min {w.min():.6f}), fluid reference not possible")
            continue
        w = np.maximum(w, -1.0)
        dat = INIS / f"t10a_w_{m}.dat"
        hdr = f"# a w  (from A0 EFT background of {m}, eft_gpresq/eft_grhoq)\n"
        dat.write_text(hdr + "".join(f"{x:.16e} {y:.16e}\n" for x, y in zip(a, w)))
        ini = harness.parse_ini(harness.PARAMS_DIR / f"{m}.ini")
        base = "hdsk_base_params.ini" if "hdsk_base" in (harness.PARAMS_DIR / f"{m}.ini").read_text() else "base_params.ini"
        extra = ""
        if m.startswith(("5_Kmouflage", "5_Kmimic")):
            extra = f"omnuh2 = {ini.get('omnuh2', '0.0')}\n"
        txt = (f"# t10a: CAMB fluid dark energy (EFTflag=0, cs2=1) with w(a) tabulated from the EFT background of {m}\n"
               f"EFTflag = 0\nuse_tabulated_w = T\nwafile = {dat}\ncs2_lam = 1\n{extra}"
               f"output_root = t10a_GR_w_{m}\nDEFAULT({harness.PARAMS_DIR / base})\n")
        (INIS / f"t10a_GR_w_{m}.ini").write_text(txt)
        print(f"{m}: wrote {dat.name} (w in [{w.min():.4f}, {w.max():.4f}]) and t10a_GR_w_{m}.ini")


if __name__ == "__main__":
    main()
