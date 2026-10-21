"""Verifier (polterddot-nonflat-gr): rerun the Kf(2) bitmask diagnostic tree (AUDIT_POLTERK) vs B0 at omk=-0.05,+0.05."""
import json, os, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness
R = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/"
P = R + "work/verify/polterK/attrib_tree"
W = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12, omk=float(sys.argv[2]), mnu=0.06, tau=0.06)
pars.InitPower.set_params(As=2.1e-9, ns=0.965)
pars.set_for_lmax(2000, lens_potential_accuracy=0)
pars.Reion.max_redshift = 50.0
cl = camb.get_results(pars).get_unlensed_scalar_cls(CMB_unit="muK", raw_cl=True)
print(json.dumps({"TT": cl[:31, 0].tolist(), "EE": cl[:31, 1].tolist(), "TE": cl[:31, 3].tolist()}))
"""
def run(tree, omk, mask=0):
    env = dict(os.environ, OMP_NUM_THREADS="1", AUDIT_POLTERK=str(mask)); env.pop("PYTHONPATH", None)
    with harness.slot(1):
        o = subprocess.run([sys.executable, "-c", W, tree, str(omk)], env=env, capture_output=True, text=True, check=True)
    return json.loads(o.stdout.strip().splitlines()[-1])
for omk in (-0.05, 0.05):
    b = run(R + "trees/B0", omk)
    for m in (0, 1, 2, 3):
        a = run(P, omk, m)
        tt = max(abs(a["TT"][l] / b["TT"][l] - 1) for l in range(2, 31))
        te = max(abs(a["TE"][l] - b["TE"][l]) / (b["TT"][l] * b["EE"][l]) ** .5 for l in range(2, 31))
        print(f"omk={omk:+.2f} mask={m} (1=octg*Kf2, 2=Edot3*Kf2)  max_l<=30 |dTT/TT|={tt:.3e} |dTE|/sqrt={te:.3e} dTT(l=2)={a['TT'][2]/b['TT'][2]-1:+.3e}")
