"""Track t4b-derivs: background-only Einstein factors (EFTeomV, EFTeomX) and their stored tau-derivative formulas,
re-evaluated offline from the get_eft_functions stencils (deriv_worker.py npz), to separate formula errors from
the integration noise of get_scale_evolution.

    python einstein_bg_offline.py NPZ_DIR [model ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deriv_analyze import dlna, load  # noqa: E402


def factors(st, a):
    Om, OmP, OmPP = st["EFTOmegaV"], st["EFTOmegaP"], st["EFTOmegaPP"]
    G4, G4P, G4PP = st["EFTGamma4V"], st["EFTGamma4P"], st["EFTGamma4PP"]
    H = st["adotoa"]
    opo = 1 + Om
    aOmP = a * OmP
    V = 0.5 / opo * (aOmP - 2 * G4 - a * G4P)
    Vdot = 0.5 * H / opo * (aOmP - 3 * a * G4P + a**2 * (OmPP - G4PP) + aOmP / opo * (-aOmP + 2 * G4 + a * G4P))
    X = 1 - G4 / opo
    Xdot = -a * H / opo * (G4P - G4 * OmP / opo)
    return V, Vdot, X, Xdot


def main():
    d_dir = Path(sys.argv[1])
    models = sys.argv[2:] or sorted(p.stem for p in d_dir.glob("*.npz"))
    print(f"{'model':28s} {'EFTeomVdot err':>15s} {'@a':>9s} {'spread':>9s} | {'EFTeomXdot err':>15s} {'@a':>9s}")
    for m in models:
        d = load(d_dir / f"{m}.npz")
        hs = d["hs"]
        x = d["xc"]
        n_h = len(hs)
        offs = np.arange(-4, 5)
        a_st = np.exp(x[None, None, :] + hs[:, None, None] * offs[None, :, None])
        st = {f: d[f"st_{f}"] for f in ("EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTGamma4V", "EFTGamma4P",
                                        "EFTGamma4PP", "adotoa")}
        V, Vdot, X, Xdot = factors(st, a_st)
        H = st["adotoa"][0, 4]
        a_on = max(d["meta"].get("pert_turn_on_results", 0.0), d["meta"].get("back_turn_on", 0.0))
        act = np.exp(x - 4 * hs.max()) > a_on
        out = []
        for F, Fd in ((V, Vdot), (X, Xdot)):
            pred = H[None, :] * dlna(F, hs)
            S = Fd[0, 4]
            sc = np.maximum.reduce([np.abs(S), np.abs(pred[1]), H * np.abs(F[0, 4]), H * 1e-3])
            err = np.abs(S - pred[1]) / sc
            spr = np.max(np.abs(pred - pred[1][None, :]), axis=0) / sc
            exc = np.where(act, np.maximum(err - 3 * spr, 0), 0)
            i = int(np.argmax(exc))
            out.append((exc[i], np.exp(x[i]), spr[i]))
        (ev, av, sv), (ex, ax, sx) = out
        print(f"{m:28s} {ev:15.2e} {av:9.2e} {sv:9.2e} | {ex:15.2e} {ax:9.2e}")


if __name__ == "__main__":
    main()
