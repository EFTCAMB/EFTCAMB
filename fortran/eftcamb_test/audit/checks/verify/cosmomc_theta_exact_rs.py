"""Verifier: is theta_MC(use_background=T) = r_s^exact(z*) / D_A, i.e. the exact-R sound horizon?"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness
os.environ["OMP_NUM_THREADS"] = "1"
TREE = sys.argv[1]
PAR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "parameters"))
model = sys.argv[2] if len(sys.argv) > 2 else "5_hdsk_TG_1"
with harness.slot(1):
    sys.path.insert(0, TREE)
    import camb
    tmp = tempfile.mkdtemp()
    out = {}
    for ub in ("F", "T"):
        ini = os.path.join(tmp, f"t{ub}.ini")
        open(ini, "w").write(f"DEFAULT({PAR}/{model}.ini)\nEFTCAMB_use_background = {ub}\nfeedback_level = 0\n")
        cp = camb.read_ini(ini, no_validate=True)
        res = camb.get_background(cp)
        ombh2, omdmh2 = cp.ombh2, cp.omch2 + cp.omnuh2
        zs = 1048 * (1 + 0.00124 * ombh2**-0.738) * (1 + (0.0783 * ombh2**-0.238 / (1 + 39.5 * ombh2**0.763))
             * (omdmh2 + ombh2) ** (0.560 / (1 + 21.1 * ombh2**1.81)))
        DA = res.angular_diameter_distance(zs) * (1 + zs)
        rs_exact = res.sound_horizon(zs)  # CAMBdata_sound_horizon (routes through EFT table when ub=T)
        out[ub] = dict(theta=res.cosmomc_theta(), rs_exact=float(rs_exact), DA=float(DA),
                       th_exact=float(rs_exact) / float(DA))
        # get_derived also
        out[ub]["rdrag"] = res.get_derived_params()["rdrag"]
for ub in "FT":
    print(ub, out[ub])
print("theta_T/theta_F-1 =", out["T"]["theta"] / out["F"]["theta"] - 1)
print("theta_T / (rs_exact_F/DA_F) -1 =", out["T"]["theta"] / out["F"]["th_exact"] - 1)
print("rs_exact T/F -1 =", out["T"]["rs_exact"] / out["F"]["rs_exact"] - 1)
print("DA T/F -1 =", out["T"]["DA"] / out["F"]["DA"] - 1)
