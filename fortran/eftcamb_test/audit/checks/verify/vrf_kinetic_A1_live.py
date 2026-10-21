"""Independent check: evaluate the timestep cache at fixed a (EFTTestStability) and compare EFT_kinetic with
36 (1+Omega)^2 A1 (gamma3=gamma4=0 identity) and with a single-h0 rebuild from the cached EFT functions."""
import ctypes, json, os, sys
from ctypes import byref, c_double, c_int
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
TREE = sys.argv[1]
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, TREE)
PAR = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
with harness.slot(1):
    import numpy as np
    import camb
    from camb.baseconfig import camblib as lib
    from camb.eftcamb import EFTCAMB_timestep_cache as C
    test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
    for m in sys.argv[2:]:
        cp = camb.read_ini(PAR + m + ".ini", no_validate=True)
        cp.EFTCAMB.EFTCAMB_skip_stability = True
        P = camb.get_background(cp).Params
        h0 = P.EFTCAMB_parameter_cache.h0_Mpc
        buf = (ctypes.c_char * (1 << 18))()
        out = []
        for a in (0.1, 0.5, 1.0):
            test(byref(c_double(a)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
            raw = bytes(buf[:8 * 400])
            g = lambda f: float(np.frombuffer(raw, np.float64, 1, getattr(C, f).offset)[0])
            O, g1, g2, g3, g4, c, H, OP = (g(f) for f in ("EFTOmegaV", "EFTGamma1V", "EFTGamma2V", "EFTGamma3V",
                                                           "EFTGamma4V", "EFTc", "adotoa", "EFTOmegaP"))
            K, A1 = g("EFT_kinetic"), g("EFTpiA1")
            Mt = 1 + O - g4
            single = 9*Mt*(4*c*Mt + 3*H**2*OP**2*a**2 + a**2*(h0**2*(3*g2**2 + 8*g1*Mt) + 6*h0*H*g2*OP))
            double = 9*Mt*(4*c*Mt + 3*H**2*OP**2*a**2 + a**2*h0**2*(3*g2**2 + 8*g1*Mt + 6*H*g2*OP))
            out.append(f"a={a}: g3={g3:.1e} g4={g4:.1e} K={K:.4e} rebuild_h0^2={double:.4e} rebuild_h0={single:.4e} "
                       f"36(1+O)^2A1={36*(1+O)**2*A1:.4e} ratio_A1/K={36*(1+O)**2*A1/K:.4f} "
                       f"ratio_A1/single={36*(1+O)**2*A1/single:.4f} cross/total={6*h0*H*g2*OP*a**2*9*Mt/single:.3f}")
        print(m); print("  " + "\n  ".join(out))
