"""Context for additional-derivs-wrong-cache: re-evaluate the ported compute_mass_factors (from
verify_addder_mass_port.py outputs) with (i) code inputs, (ii) analytic derivatives, (iii) code inputs with, and (iv) analytic derivatives with,
EFTcPP replaced by the dimensionally consistent d(EFTcP)/da = [cdotdot + (3H - Hdot/H) cdot + 2 H^2 c]/(a^2 H^2).
Usage: python verify_addder_cpp_context.py PORT_OUTDIR
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_addder_mass_port as vp  # noqa: E402

outdir = Path(sys.argv[1])
src, _ = vp.fortran_mass_factors()
old = "EFTcPP = (c['EFTcdotdot']/c['adotoa'] + (3. - c['Hdot']/c['adotoa'])*c['EFTcdot'] + 2.*c['adotoa']*c['EFTc'])/(a**2*c['adotoa']**2)"
assert old in src, "EFTcPP statement not found"
new = ("EFTcPP = (c['EFTcdotdot'] + (3.*c['adotoa'] - c['Hdot']/c['adotoa'])*c['EFTcdot'] + "
       "2.*c['adotoa']**2*c['EFTc'])/(a**2*c['adotoa']**2)")
src_fixed = src.replace(old, new)
for f in sorted(outdir.glob("A0_*.npz")):
    d = np.load(f)
    a, ton, rate, h0 = d["a"], float(d["ton"]), float(d["rate"]), float(d["h0_Mpc"])
    if a.min() > ton:
        pass
    used, ana = vp.todict(d["used"]), vp.todict(d["ana"])
    H = used["adotoa"]
    thr = -(rate * H / a) ** 2
    act = a > ton
    if not act.any():
        print(f"{f.stem[3:]:24s} mass check never active (turn-on {ton})")
        continue
    corr = dict(used)
    for k in vp.REPLACED:
        corr[k] = ana[k]
    res = []
    for s, c in ((src, used), (src, corr), (src_fixed, used), (src_fixed, corr)):
        m1, m2 = vp.mass_factors(s, c, a, h0)
        uns = bool(np.any(((m1 < thr) | (m2 < thr)) & act))
        res.append(f"{'UNSTABLE' if uns else 'stable  '} min mu1/H2={np.nanmin((m1 * (a / H) ** 2)[act]):+.3e}")
    print(f"{f.stem[3:]:24s} code-inputs: {res[0]} | analytic derivs: {res[1]} | EFTcPP fixed, code inputs: {res[2]} | EFTcPP fixed + analytic: {res[3]}")
