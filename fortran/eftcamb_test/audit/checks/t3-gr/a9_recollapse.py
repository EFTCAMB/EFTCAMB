"""t3-gr (a): EFTflag=0 dtauda restructure: sqrt(3/grhoa2) is evaluated before the grhoa2<=0 (recollapse) check.

Closed model with Omega_m=0.05, Omega_k=-1.05, Omega_L=2: grhoa2(a=0.5) < 0. B0/A0 (release) should raise the CAMB
error 'Universe stops expanding before today'; AD (FPE traps) shows whether the invalid sqrt traps first.
"""
from common import harness, run

ov = {"hubble": "70", "ombh2": "0.01", "omch2": "0.0145", "omk": "-1.05", "omnuh2": "0", "massive_neutrinos": "0",
      "massless_neutrinos": "3.046", "nu_mass_eigenstates": "0", "num_redshiftwindows": "0", "get_tensor_cls": "F",
      "do_lensing": "F", "get_transfer": "F", "l_max_scalar": "200"}
for t in ("B0", "A0", "AD"):
    p = run(t, f"{t}_recollapse", overrides=ov, sub="a9", timeout=600)
    st = harness.json.loads(p.with_suffix(".json").read_text())
    tail = [ln for ln in st.get("log_tail", []) if "rror" in ln or "FPE" in ln or "Floating" in ln or "stops" in ln]
    print(t, st["status"], tail[-4:])
