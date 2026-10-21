"""Verifier: evolve_metric_h T vs F on TG_1 variants (cubic off; massless nu) + scale evolution at a few k."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "evolve_metric_h"
OUT.mkdir(parents=True, exist_ok=True)
INI = harness.PARAMS_DIR / "5_hdsk_TG_1.ini"
MLESS = {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0"}
VAR = {
    "nocubic": {"Horndeskic30": "0"},
    "masslessnu": MLESS,
}
EK = [1e-3, 0.05, 1.0, 30.0]
jobs = []
for v, ov in VAR.items():
    for h in ("F", "T"):
        jobs.append((f"{v}_{h}", {**ov, "EFTCAMB_evolve_metric_h": h}, None))
for h in ("F", "T"):
    jobs.append((f"evol_{h}", {"EFTCAMB_evolve_metric_h": h}, EK))


def run(j):
    name, ov, ek = j
    return name, harness.run_model("A0", INI, OUT / f"{name}.npz", overrides=ov, evol_k=ek, threads=1)


with ThreadPoolExecutor(3) as ex:
    st = dict(ex.map(run, jobs))
for n, s in st.items():
    print(n, s["status"])

ells = [2, 5, 10, 30, 100]
ks = [0.1, 1.0, 10.0, 50.0, 100.0]


def cmp(a, b):
    ra, rb = harness.load(OUT / f"{a}.npz"), harness.load(OUT / f"{b}.npz")
    n = min(ra["cl_unlensed_scalar"].shape[0], rb["cl_unlensed_scalar"].shape[0])
    tt = rb["cl_unlensed_scalar"][:n, 0] / ra["cl_unlensed_scalar"][:n, 0] - 1
    km = min(ra["pk_k"].max(), rb["pk_k"].max())
    pa = np.interp(np.log(ks), np.log(ra["pk_k"]), np.log(ra["pk"][0]))
    pb = np.interp(np.log(ks), np.log(rb["pk_k"]), np.log(rb["pk"][0]))
    pk = np.where(np.array(ks) <= km, np.exp(pb - pa) - 1, np.nan)
    print(f"{b} vs {a}: dTT " + " ".join(f"l{e}:{tt[e]:+.2e}" for e in ells) + " | dP " +
          " ".join(f"k{k:g}:{x:+.2e}" for k, x in zip(ks, pk)))


for v in VAR:
    if st[f"{v}_F"]["status"] == "ok" and st[f"{v}_T"]["status"] == "ok":
        cmp(f"{v}_F", f"{v}_T")
if st["evol_F"]["status"] == "ok" and st["evol_T"]["status"] == "ok":
    cmp("evol_F", "evol_T")
    F, T = harness.load(OUT / "evol_F.npz"), harness.load(OUT / "evol_T.npz")
    print("evol fields:", [k for k in F if k.startswith("evol_")])
    a = F["evol_a"]
    for fld in ("evol_z", "evol_sigma", "evol_etak", "evol_delta_phi", "evol_clxc"):
        if fld not in F:
            continue
        for i, k in enumerate(EK):
            f, t = F[fld][i], T[fld][i]
            for ai in (1e-3, 1e-2, 0.1, 0.5, 1.0):
                j = np.argmin(abs(np.log(a / ai)))
                print(f"{fld} k={k:g} a={ai:g}: F={f[j]:+.4e} T={t[j]:+.4e} rel={(t[j]-f[j])/(abs(f[j])+1e-300):+.2e}")
