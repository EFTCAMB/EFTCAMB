"""Probe how cosmomc_theta of 5_hdsk_TG_1 responds to the input H0 (use_background F/T), to see whether the
theta -> H0 mapping (set_H0_for_theta) is usable for this model. Runs inside a harness CPU slot."""
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import harness  # noqa: E402

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
OUT = os.path.join(ROOT, "work", "verify", "cosmomc_theta_usebg")
model = sys.argv[1] if len(sys.argv) > 1 else "5_hdsk_TG_1"
extra = sys.argv[2:]  # extra ini lines key=value
sys.path.insert(0, os.path.join(ROOT, "trees", "A0"))
os.nice(10)
with harness.slot(1):
    import camb
    for ub in ("F", "T"):
        ini = os.path.join(OUT, f"probe_{model}_ub{ub}.ini")
        with open(ini, "w") as f:
            f.write(f"DEFAULT(/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/{model}.ini)\n"
                    f"EFTCAMB_use_background = {ub}\nfeedback_level = 0\n" + "".join(e + "\n" for e in extra))
        for H0 in (65.0, 68.0, 70.86711257, 74.0, 78.0):
            cp = camb.read_ini(ini, no_validate=True)
            cp.H0 = H0
            try:
                res = camb.get_background(cp)
                print(f"ub={ub} H0_in={H0:8.4f} theta={res.cosmomc_theta():.9f} H(z=0)={float(res.hubble_parameter(0.0)):.5f}"
                      f" h0_cache={res.Params.eft_par_cache.h0 if hasattr(res.Params, 'eft_par_cache') else 'na'}")
            except Exception as e:  # noqa: BLE001
                print(f"ub={ub} H0_in={H0:8.4f} ERROR {e!r}")
