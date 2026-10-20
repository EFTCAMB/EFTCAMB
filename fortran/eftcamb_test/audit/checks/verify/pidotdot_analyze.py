"""Analyse dumps from pidotdot_dump.py: derivative-consistency of EFTISW/PsiDot/sigmadotdot (affected by pi'' bug)
against PhiDot (control), A0 vs FX, plus direct check that FX-A0 equals the analytically missing pi'' terms."""
import sys

import numpy as np

D = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
     "scratchpad/audit/work/verify/pidotdot/dbg/")
sl = slice(10, -10)


def resid(val, ref):
    return np.max(np.abs(val[sl] - ref[sl])) / np.max(np.abs(ref[sl]))


for m in sys.argv[1:]:
    dd = {t: np.load(f"{D}{m}_{t}.npz") for t in ("A0", "FX")}
    h0 = float(dd["A0"]["h0_mpc"])
    print(f"== {m}")
    for ik, k in enumerate(dd["A0"]["ks"]):
        line = []
        for t, d in dd.items():
            c = {f: d[f][ik] for f in d.files if d[f].ndim == 2}
            tau = c["tau"]
            g = lambda f: np.gradient(c[f], tau, edge_order=2)  # noqa: E731
            hd = resid(c["Hdot"], g("adotoa"))
            line.append(f"{t}: ISW {resid(c['EFTISW'], g('EFTLensing')):.1e} PsiDot {resid(c['PsiDot'], g('Psi')):.1e} "
                        f"sigdd {resid(c['sigmadotdot'], g('sigmadot')):.1e} PhiDot(ctrl) {resid(c['PhiDot'], g('Phi')):.1e} "
                        f"Hdot {hd:.0e}")
        print(f"  k={k:<6} " + " | ".join(line))
        # direct: FX - A0 vs missing terms computed from A0 cache
        a, f = ({x: dd[t][x][ik] for x in dd[t].files if dd[t][x].ndim == 2} for t in ("A0", "FX"))
        same_traj = np.max(np.abs(a["pidotdot"] - f["pidotdot"])) / np.max(np.abs(a["pidotdot"]) + 1e-300)
        miss_nd = a["k"] * (a["EFTGamma4V"] + 2 * a["EFTGamma5V"]) * a["pidotdot"] / (1 + a["EFTOmegaV"])
        miss_isw = miss_nd / (a["EFTeomX"] * h0)
        d_isw = f["EFTISW"] - a["EFTISW"]
        d_nd = f["EFTeomNdot"] - a["EFTeomNdot"]
        rel = np.max(np.abs(d_isw - miss_isw)) / max(np.max(np.abs(miss_isw)), 1e-300)
        reln = np.max(np.abs(d_nd - miss_nd)) / max(np.max(np.abs(miss_nd)), 1e-300)
        print(f"           pi'' A0 vs FX rel {same_traj:.1e}; max|missing ISW|/max|ISW| = "
              f"{np.max(np.abs(miss_isw))/np.max(np.abs(a['EFTISW'])):.2e}; "
              f"(FX-A0) vs missing: ISW {rel:.1e}, Ndot {reln:.1e}; "
              f"max|dM|/max|M| = {np.max(np.abs(f['EFTeomM']-a['EFTeomM']))/np.max(np.abs(a['EFTeomM'])):.1e}")
