"""Independent check: an EFT model loaded via camb.read_ini, then set_params(cp=pars, ns=...) (no cosmology keys,
so set_cosmology is not called). Does the EFT model survive? Also direct initialize_parameters(cp, {}) call."""
import json, os, subprocess, sys
REPO = "/Users/marco/IW/eftcamb_update/EFTCAMB"
TREES = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/"
WORKER = r"""
import sys, json, os, faulthandler; faulthandler.enable()
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[3] + "/fortran/eftcamb_test/audit/checks")
import harness
with harness.slot(1):
    import camb
    out = {}
    ini = sys.argv[3] + "/fortran/eftcamb_test/parameters/" + sys.argv[2]
    p = camb.read_ini(ini)
    if p.EFTCAMB.EFTflag == 0:  # B1 read_ini ignores EFT keys: init from flattened ini like harness
        flat = harness.parse_ini(ini)
        p.EFTCAMB.initialize_parameters(p, flat, print_header=False)
    p.set_for_lmax(500); p.WantTransfer = False
    out["flag_before"] = int(p.EFTCAMB.EFTflag)
    out["pars_before"] = {k: v for k, v in p.EFTCAMB.read_parameters().items() if isinstance(v, (int, float))}
    out["s8_before"] = float(camb.get_results(p).get_cmb_power_spectra(CMB_unit="muK")["total"][2,0])
    ns0 = p.InitPower.ns
    p2 = camb.set_params(cp=p, ns=ns0)  # same ns: physics should be identical
    out["flag_after"] = int(p2.EFTCAMB.EFTflag)
    out["same_object"] = p2 is p
    out["s8_after"] = float(camb.get_results(p2).get_cmb_power_spectra(CMB_unit="muK")["total"][2,0])
    out["H0_after"] = float(p2.H0)
print("JSON" + json.dumps(out))
"""
def run(tree, ini):
    env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
    r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", WORKER, TREES + tree, ini, REPO], capture_output=True, text=True, env=env)
    print("RC", r.returncode, r.stderr[-3000:])
    l = [x for x in r.stdout.splitlines() if x.startswith("JSON")]
    if not l: return {"error": (r.stdout + r.stderr)[-6000:]}
    return json.loads(l[0][4:])
for tree in sys.argv[2:] or ("A0", "B1"):
    for ini in ("4_Des_F_R_LCDM_1.ini", "2_PEFT_Omega_const_1.ini"):
        o = run(tree, ini)
        if "error" in o: print(tree, ini, o["error"]); continue
        print(f"{tree} {ini}: EFTflag {o['flag_before']} -> {o['flag_after']}  TT(l=2) {o['s8_before']:.6f} -> {o['s8_after']:.6f}  H0={o['H0_after']}")
