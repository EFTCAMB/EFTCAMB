"""Verifier: halofit w0wa for EFT. Independent attribution checks.
usage: python v_nl_w0wa.py TREE w_value version
Compares P_nl/P_lin of pure-EFT designer wCDM (EFTw0=w) vs CAMB fluid (w, cs2=1), and prints effective w0wa."""
import os, sys, tempfile
from pathlib import Path
os.environ["OMP_NUM_THREADS"] = "2"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa
import numpy as np
TREE, W, VER = sys.argv[1], float(sys.argv[2]), sys.argv[3]
sys.path.insert(0, TREE)
import camb
from camb import model
P = Path(__file__).resolve().parents[3] / "parameters"
tmp = Path(tempfile.mkdtemp())
COMMON = "get_tensor_cls = F\nget_scalar_cls = F\nnum_redshiftwindows = 0\nfeedback_level = 0\nEFTCAMB_write_background = F\n"

def run(name, extra):
    ini = tmp / f"{name}.ini"
    ini.write_text(f"DEFAULT({P / name}.ini)\n{COMMON}{extra}")
    cp = camb.read_ini(str(ini))
    keys = harness.parse_ini(str(ini))
    if cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
    cp.set_matter_power(redshifts=[1.0, 0.0], kmax=10.0, nonlinear=True)
    cp.NonLinear = model.NonLinear_pk
    cp.NonLinearModel.set_params(halofit_version=VER)
    r = camb.get_results(cp)
    k, z, pl = r.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    _, _, pn = r.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False)
    w0wa = None
    try:
        w0wa = cp.EFTCAMB.get_effective_w0wa()
    except Exception as e:
        w0wa = repr(e)
    return cp, k, z, pl, pn / pl, w0wa

with harness.slot(2):
    cpe, k, z, ple, re_, we = run("2_PEFT_wDE_wCDM_1", f"EFTw0 = {W}\n")
    cpg, k2, _, plg, rg, wg = run("1_EFT_GR", f"w = {W}\ncs2_lam = 1\n")
sel = (k > 0.05) & (k < 5)
print(f"tree={Path(TREE).name} w={W} ver={VER} EFTflag_e={cpe.EFTCAMB.EFTflag} EFTflag_g={cpg.EFTCAMB.EFTflag} "
      f"w0wa_eft={we} fluid_w={cpg.DarkEnergy.w}")
for i, zz in enumerate(z):
    dl = np.max(np.abs(ple[i, sel] / plg[i, sel] - 1)); dr = re_[i, sel] / rg[i, sel] - 1; j = np.argmax(np.abs(dr))
    print(f"  z={zz:.1f} lin={dl:.2e} ratio={dr[j]:+.3e} at k={k[sel][j]:.2f}")
