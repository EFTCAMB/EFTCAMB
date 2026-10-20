"""Independent check for stability-gpinu-tot-tripled: finite-difference dHdot/dtau and dHdotdot/dtau from the
python bridge (single-reset loop) vs the analytic Hdotdot of (a) the bridge and (b) the stability routine
(EFTTestStability) in A0 and fixnu, at the same a. Also direct ThermalNu-free test: stability gpinu_tot vs
3*gpresm-derived massive-nu pressure (gpresm_t - (grhor+grhog)/3)."""
import ctypes, os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness
tree = sys.argv[1]
sys.path.insert(0, tree)
import numpy as np
with harness.slot(1):
    import camb
    from camb.baseconfig import camblib as lib
    from camb.eftcamb import EFTCAMB_timestep_cache
    from ctypes import byref, c_double, c_int
    ini = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/" + sys.argv[2] + ".ini"
    cp = camb.read_ini(ini, no_validate=True)
    cp.EFTCAMB.EFTCAMB_skip_stability = True
    res = camb.get_background(cp)
    P = res.Params
    a0 = np.array([1e-6, 1e-4, 1e-3, 1e-2, 1e-1])
    h = 1e-4
    grid = np.concatenate([a0 * np.exp(-h), a0, a0 * np.exp(h)])
    _, br = cp.EFTCAMB.get_eft_functions(res, grid)
    n = len(a0)
    def sl(f, i): return np.asarray(br[f])[i * n:(i + 1) * n]
    H = sl("adotoa", 1)
    # d/dtau = a*H d/da = H d/dln a
    fd_Hdd = H * (sl("Hdot", 2) - sl("Hdot", 0)) / (2 * h)
    fd_Hddd = H * (sl("Hdotdot", 2) - sl("Hdotdot", 0)) / (2 * h)
    test = getattr(lib, "__eftcamb_stability_MOD_eftteststability"); test.restype = c_int
    buf = (ctypes.c_char * (1 << 18))()
    names = ["gpinu_tot", "gpresm_t", "grhor_t", "grhog_t", "Hdotdot", "Hdotdotdot"]
    off = {k: getattr(EFTCAMB_timestep_cache, k).offset for k in names}
    st = {k: [] for k in names}
    for x in a0:
        test(byref(c_double(x)), byref(c_double(10.0)), byref(P.EFTCAMB), byref(P.EFTCAMB_parameter_cache), buf)
        raw = bytes(buf[: max(off.values()) + 8])
        for k in names:
            st[k].append(float(np.frombuffer(raw, dtype=np.float64, count=1, offset=off[k])[0]))
    st = {k: np.array(v) for k, v in st.items()}
    nu_from_gpres = st["gpresm_t"] - (st["grhor_t"] + st["grhog_t"]) / 3
print("tree", tree, "model", sys.argv[2])
print(" a       stab_gpinu/gpres_nu  (brHdd-FD)/FD  (stHdd-FD)/FD  (brHddd-FD)/FD  (stHddd-FD)/FD")
for i, x in enumerate(a0):
    print(f"{x:8.1e} {st['gpinu_tot'][i]/nu_from_gpres[i]:10.5f} {(sl('Hdotdot',1)[i]-fd_Hdd[i])/abs(fd_Hdd[i]):14.2e}"
          f" {(st['Hdotdot'][i]-fd_Hdd[i])/abs(fd_Hdd[i]):14.2e} {(sl('Hdotdotdot',1)[i]-fd_Hddd[i])/abs(fd_Hddd[i]):14.2e}"
          f" {(st['Hdotdotdot'][i]-fd_Hddd[i])/abs(fd_Hddd[i]):14.2e}")
