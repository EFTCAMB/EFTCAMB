"""Track 5: focused side runs, one CAMB process at a time (restartable: skips existing outputs).

  fr/*      f(R) designer B0=1e-4 (4_Des_F_R_LCDM_2): GR-threshold convergence (1e-10, 1e-11, skip RGR)
  hdsk/*    Horndeski TG: evolve_delta_phi / evolve_metric_h / use_background toggles
  pidfix/*  Horava models on the tree with the EFTpiDfunction fix (D1 + k^2 D2)
Usage: python queue_side.py [job,...]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "t5-robust"
P = harness.PARAMS_DIR
FR2 = P / "4_Des_F_R_LCDM_2.ini"
FR1 = P / "4_Des_F_R_LCDM_1.ini"
PIDFIX = W / "A0_pidfix"
PIDD = W / "A0_pidd"
JOBS = [
    # f(R) B0=1e-4: pi(a) at late times for the default turn-on (RGR a=0.093) vs an early turn-on (threshold 1e-12)
    ("frevol/def", "A0", FR2, {}, {}, [0.1, 1.0, 2.0, 5.0]),
    ("frevol/grt1e-12", "A0", FR2, {"EFTCAMB_GR_threshold": "1e-12"}, {}, [0.1, 1.0, 2.0, 5.0]),
    ("frevol2/pidfix", PIDFIX, FR2, {}, {}, [0.1, 1.0, 2.0, 5.0]),
    ("frevol2/AB2", "A0", FR2, {"accuracy_boost": "2"}, {}, [0.1, 1.0, 2.0, 5.0]),
    ("frevol2/pto0.2", "A0", FR2, {"EFTCAMB_turn_on_time": "0.2"}, {}, [0.1, 1.0, 2.0, 5.0]),
    # stale pi'' read in the Einstein factors (EFT_T5_PIDD: 0 stale = release, 1 zero, 2 previous derivs value)
    *[(f"pidd/TG2_mode{m}", PIDD, P / "5_hdsk_TG_2.ini", {}, {"EFT_T5_PIDD": m}) for m in ("0", "1", "2")],
    # Horndeski toggles (TG_1 = delta_phi T, TG_2 = delta_phi F, TG_3 = other parameters, delta_phi T)
    ("hdsk/TG3_dphiF", "A0", P / "5_hdsk_TG_3.ini", {"EFTCAMB_evolve_delta_phi": "F"}, {}),
    ("hdsk/TG1_ubT", "A0", P / "5_hdsk_TG_1.ini", {"EFTCAMB_use_background": "T"}, {}),
    ("hdsk/TG2_ubT", "A0", P / "5_hdsk_TG_2.ini", {"EFTCAMB_use_background": "T"}, {}),
    ("hdsk/TG1_mhT", "A0", P / "5_hdsk_TG_1.ini", {"EFTCAMB_evolve_metric_h": "T"}, {}),
    ("pidd/TG3_dphiF_mode2", PIDD, P / "5_hdsk_TG_3.ini", {"EFTCAMB_evolve_delta_phi": "F"}, {"EFT_T5_PIDD": "2"}),
    # IC types 2/3 with a turn-on deep in radiation domination (where the pi series is meant to apply)
    ("icx/Omc2_ic1_pto1e-8", "A0", P / "2_PEFT_Omega_const_2.ini", {"EFT_IC_type": "1", "EFTCAMB_turn_on_time": "1e-8"}, {}),
    ("icx/Omc2_ic2_pto1e-8", "A0", P / "2_PEFT_Omega_const_2.ini", {"EFT_IC_type": "2", "EFTCAMB_turn_on_time": "1e-8"}, {}),
    ("icx/Omc2_ic3_pto1e-8", "A0", P / "2_PEFT_Omega_const_2.ini", {"EFT_IC_type": "3", "EFTCAMB_turn_on_time": "1e-8"}, {}),
    # Horava with the pi' initial-condition fix
    *[(f"pidfix_h/{m}", PIDFIX, P / f"{m}.ini", {}, {}) for m in (
        "5_Horava_eta_1", "5_Horava_eta_3", "5_Horava_SolSyst_1", "5_Horava_eta_2", "5_Horava_eta_4",
        "5_Horava_SolSyst_2", "5_Horava_SolSyst_3")],
    ("fr/LCDM1_grt1e-12", "A0", FR1, {"EFTCAMB_GR_threshold": "1e-12", "EFTCAMB_turn_on_time": "1e-3"}, {}),
    ("pidd/frLCDM2_mode2", PIDD, FR2, {}, {"EFT_T5_PIDD": "2"}),
]


def main() -> None:
    only = sys.argv[1].split(",") if len(sys.argv) > 1 else None
    for job in JOBS:
        name, tree, ini, ov, env = job[:5]
        ks = job[5] if len(job) > 5 else None
        if only and not any(name.startswith(o) for o in only):
            continue
        out = W / f"{name}.npz"
        if out.with_suffix(".json").exists():
            continue
        st = harness.run_model(tree, ini, out, overrides=ov, threads=2, timeout=2400, extra_env=env, evol_k=ks)
        print(name, st["status"], st["elapsed"], flush=True)
    print("SIDE DONE", flush=True)


if __name__ == "__main__":
    main()
