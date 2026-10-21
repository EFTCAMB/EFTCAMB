"""Verifier: term-by-term decomposition of EFT_kinetic for 3_RPH_MB_OmegaDE_1 from the code's own cache fields,
compared with the cancellation-free form 18 M^3 H^2 (alphaK + 6 alphaB^2), M = 1 + OmegaV."""
import os, sys, json, ctypes
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
TREE = sys.argv[1] if len(sys.argv) > 1 else str(harness.TREES["A0"])
sys.path.insert(0, TREE)
from ctypes import byref, c_double, c_int
import numpy as np
with harness.slot(1):
    import camb
    from camb.baseconfig import camblib as lib
    from camb.eftcamb import EFTCAMB_timestep_cache as C
    cp = camb.read_ini(str(harness.PARAMS_DIR / "3_RPH_MB_OmegaDE_1.ini"), no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    P = camb.get_background(cp).Params
    h0 = P.EFTCAMB_parameter_cache.h0_Mpc
    test = lib.__eftcamb_stability_MOD_eftteststability; test.restype = c_int
    buf = (ctypes.c_char * (1 << 18))()
    print("h0_Mpc", h0)
    for a in (1e-9, 1e-7, 3.2e-7, 1e-6, 6e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0):
        test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
        raw = bytes(buf[:8 * 300])
        g = lambda f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0])
        OV, OP, c, g1, g2, g4, H = g("EFTOmegaV"), g("EFTOmegaP"), g("EFTc"), g("EFTGamma1V"), g("EFTGamma2V"), g("EFTGamma4V"), g("adotoa")
        M = 1 + OV - g4
        t_c = 4 * c * M
        t_g1 = 8 * g1 * M * a**2 * h0**2
        t_rest = 3 * H**2 * OP**2 * a**2 + a**2 * h0**2 * 3 * g2**2 + 6 * a**2 * h0 * H * g2 * OP
        code_py = 9 * M * (t_c + 3 * H**2 * OP**2 * a**2 + a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * M) + 6 * H * g2 * OP))
        aK, aB = g("alphaK"), g("alphaB")
        ana = 18 * M**3 * H**2 * (aK + 6 * aB**2)
        ana_split = 9 * M * (2 * aK * M**2 * H**2 + 12 * aB**2 * M**2 * H**2)
        print(f"a={a:.1e} kin_code={g('EFT_kinetic'): .4e} code_py={code_py: .4e} ana={ana: .4e} "
              f"|4cM|={abs(t_c):.3e} (4cM+8g1..)={t_c + t_g1: .3e} vs 2aK M^2H^2={2*aK*M**2*H**2: .3e} "
              f"rest={t_rest: .3e} vs 12aB^2M^2H^2={12*aB**2*M**2*H**2: .3e} aB={aB: .3e} aK={aK: .3e} OV={OV:.3e} Meff2={g('Meff2'):.4e}")
