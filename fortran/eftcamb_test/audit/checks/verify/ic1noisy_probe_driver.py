"""Driver for ic1noisy_probe.py: runs probe cases one at a time in a subprocess holding a harness CPU slot."""

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "probe"
OUT.mkdir(parents=True, exist_ok=True)
CASES = {
    "km_A0_1e-8": ("A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-8"}, 1e-8),
    "km_A0_1.01e-8": ("A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1.01e-8"}, 1.01e-8),
    "km_A0_0.97e-8": ("A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "0.97e-8"}, 0.97e-8),
    "km_A0_1e-7": ("A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-7"}, 1e-7),
    "km_A0_1e-9": ("A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-9"}, 1e-9),
    "km_A0_1e-8_np5000": ("A0", "5_Kmimic_1", {"EFTCAMB_turn_on_time": "1e-8",
                                               "model_background_num_points": "5000"}, 1e-8),
    "tg2_A0_1e-8": ("A0", "5_hdsk_TG_2", {}, 1e-8),
    "tg2_A0_2e-8": ("A0", "5_hdsk_TG_2", {"EFTCAMB_turn_on_time": "2e-8"}, 2e-8),
}
for name in sys.argv[1:] or list(CASES):
    tree, model, ov, a_on = CASES[name]
    out = OUT / f"{name}.json"
    if out.exists():
        continue
    wrapper = harness._wrapper_ini(harness.PARAMS_DIR / f"{model}.ini", OUT, name, ov)
    half = harness._wrapper_ini(harness.PARAMS_DIR / f"{model}.ini", OUT, name + "_half",
                                {**ov, "EFTCAMB_turn_on_time": repr(0.5 * a_on)})
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1), open(OUT / f"{name}.log", "w") as f:
        rc = subprocess.run([sys.executable, str(HERE / "ic1noisy_probe.py"), str(harness.tree_path(tree)),
                             str(half), str(wrapper), repr(a_on), str(out)], cwd=OUT, env=env, stdout=f,
                            stderr=subprocess.STDOUT, timeout=500, preexec_fn=harness._nice).returncode
    print(name, "rc", rc)
