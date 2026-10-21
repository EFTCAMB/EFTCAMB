"""t3-gr (b): Horndeski (EFTflag=5, model 7: G2 = X - c0(phi) + c1 X^2, G4 = 1/2 + c3 + c4 X) at its LCDM point:
c0 = constant (shot to Omega_DE), c3 = c1 = c4 = 0, phi_dot_ini = 0 -> a cosmological constant. Compared with
EFTflag=0 at the same cosmology (hdsk_base_params.ini). l_max reduced from the hdsk base (10000) to keep it cheap.
"""

from concurrent.futures import ThreadPoolExecutor

from common import INIS, PARAMS, RUNS, fmt, harness

SUB = RUNS / "b9"
LOW = {"l_max_scalar": "2500", "k_eta_max_scalar": "5000", "l_max_tensor": "1500", "k_eta_max_tensor": "3000"}
HDSK = PARAMS / "5_hdsk_TG_2.ini"
# c0 = 8 pi G rho_Lambda in Mpc^-2 (G4 = 1/2 means 8 pi G = 1): grhov = 1.16537e-07 for this cosmology
GR = {"Horndeski_freefunc0_model": "1", "Horndeskic00": "1.16537e-07", "Horndeski_freefunc3_model": "0",
      "Horndeski_shooting": "F", "EFTCAMB_skip_RGR": "F", "EFTCAMB_skip_stability": "F"}
JOBS = {
    "B0_ref": ("B0", INIS / "t10a_GR_hdsk.ini", LOW, {}),
    "A0g_ref": ("A0", INIS / "t10a_GR_hdsk.ini", LOW, {"Accuracy.IntTolBoost": 2.0}),
    "hdsk_GR": ("A0", HDSK, {**LOW, **GR}, {}),
    "hdsk_GR_forced": ("A0", HDSK, {**LOW, **GR, "EFTCAMB_skip_RGR": "T", "EFTCAMB_skip_stability": "T"}, {}),
    "hdsk_GR_shoot": ("A0", HDSK, {**LOW, **GR, "Horndeski_shooting": "T"}, {}),
}


def one(item):
    n, (tree, ini, ov, ps) = item
    out = SUB / f"{n}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model(tree, ini, out, overrides=ov, python_set=ps, threads=1, timeout=540)
    return n, out


with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(one, JOBS.items()))
for n in ("hdsk_GR", "hdsk_GR_forced", "hdsk_GR_shoot"):
    st = harness.json.loads(res[n].with_suffix(".json").read_text())
    if st["status"] != "ok":
        print(n, st["status"], [ln for ln in st.get("log_tail", []) if "rror" in ln or "nstab" in ln or "fail" in ln][-3:])
        continue
    for r in ("A0g_ref", "B0_ref"):
        print(f"{n:15s} vs {r:8s}", fmt(harness.compare(res[n], res[r])))
