"""(iii) Horndeski test models (5_hdsk_TG_*, shipped with EFTCAMB_skip_stability = T): verdict with the stability check
forced on, as a function of EFTCAMB_stability_time (default 1e-10; perturbation turn-on of these inis = 1e-8)."""
import sys

sys.path.insert(0, ".")
import t6lib  # noqa: E402
import harness  # noqa: E402

for m in sys.argv[1:] or ["5_hdsk_TG_1", "5_hdsk_TG_3"]:
    row = []
    for st in ("1e-10", "2.5e-10", "1e-9", "1e-8"):
        r = t6lib.probe("A0", harness.PARAMS_DIR / f"{m}.ini",
                        overrides={"EFTCAMB_skip_stability": "F", "EFTCAMB_stability_time": st, "feedback_level": "0"},
                        timeout=600)
        row.append(f"stability_time={st}: {'stable' if r['bg_error'] is None else 'UNSTABLE'} ({r['status']})")
    print(m, " | ".join(row))
