"""Independent check of hdsk-h0-not-propagated from stored reference npz (A0 and B1): H0 from EFT cache adotoa(a=1),
bg_H(z=0), and sigma8 recomputed by integrating stored linear P(k) [1/Mpc units] with R = 8/h for both h values."""
import json, sys
import numpy as np
from scipy.integrate import simpson
R = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/runs"
c = 299792.458
def sig(k, pk, Rm):
    x = k * Rm
    W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    lk = np.log(k)
    return np.sqrt(simpson(k**3 * pk * W**2 / (2 * np.pi**2), x=lk))
for tree in ("A0", "B1"):
    for m in ("5_hdsk_TG_1", "5_hdsk_TG_2", "5_hdsk_TG_3"):
        d = np.load(f"{R}/{tree}/{m}.npz", allow_pickle=True)
        meta = json.loads(str(d["meta"]))
        keys = [k for k in d.files if "adotoa" in k.lower()]
        H0_eft = None
        if "eft_adotoa" in d.files:
            a = d["eft_a_grid"]; H0_eft = d["eft_adotoa"][-1] / a[-1] * c
        iz = int(np.argmin(np.abs(d["bg_z"])))
        bgH = d["bg_H"][iz]
        H0_in = meta.get("eft_params", {}).get("H0") if isinstance(meta.get("eft_params"), dict) else None
        k = d["pk_k"]; zi = int(np.argmin(np.abs(d["pk_z"]))); pk = d["pk"][zi]
        s_in = sig(k, pk, 8 / 0.7086711257); s_true = sig(k, pk, 8 / (bgH / 100))
        print(f"{tree} {m}: bg_z[{iz}]={d['bg_z'][iz]:.3g} bg_H={bgH:.4f} H0_from_eft_adotoa={H0_eft} "
              f"sigma8_stored={np.ravel(d['sigma8'])[-1]:.4f} sig(8/h_in)={s_in:.4f} sig(8/h_true)={s_true:.4f} "
              f"pk_z={d['pk_z'][zi]} kmax={k[-1]:.2f}")
