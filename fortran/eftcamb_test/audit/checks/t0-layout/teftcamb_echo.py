"""Track 0: python -> Fortran echo of TEFTCAMB fields. Set distinctive values through the ctypes fields, let
Fortran print them (EFTCAMB.feedback -> print_EFTCAMB_flags) and parse the printout.

Usage: python teftcamb_echo.py [--tree PATH]     (spawns itself in a subprocess to capture Fortran stdout)
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
VALUES = {
    "EFTCAMB_feedback_level": 7, "EFTCAMB_back_turn_on": 3.25e-9, "EFTCAMB_pert_turn_on": 0.0123,
    "EFTCAMB_GR_threshold": 4.5e-7, "EFTCAMB_stability_time": 6.5e-11, "EFTCAMB_stability_threshold": 0.375,
    "EFTCAMB_effective_w0wa": True, "EFT_ghost_math_stability": True, "EFT_mass_math_stability": False,
    "EFT_ghost_stability": False, "EFT_gradient_stability": True, "EFT_mass_stability": True,
    "EFT_mass_stability_rate": 12.5, "EFT_additional_priors": False, "EFT_positivity_bounds": False,
    "EFT_minkowski_limit": True,
}
PATTERNS = {
    "EFTCAMB_feedback_level": r"EFTCAMB_feedback_level\s*=\s*(\S+)",
    "EFTCAMB_back_turn_on": r"EFTCAMB_back_turn_on\s*=\s*(\S+)",
    "EFTCAMB_pert_turn_on": r"EFTCAMB_pert_turn_on\s*=\s*(\S+)",
    "EFTCAMB_GR_threshold": r"EFTCAMB_GR_threshold\s*=\s*(\S+)",
    "EFTCAMB_stability_time": r"EFTCAMB_stability_time\s*=\s*(\S+)",
    "EFTCAMB_stability_threshold": r"EFTCAMB_stability_threshold\s*=\s*(\S+)",
    "EFTCAMB_effective_w0wa": r"EFTCAMB_effective_w0wa\s*=\s*(\S+)",
    "EFT_ghost_math_stability": r"Math ghost stability\s*=\s*(\S+)",
    "EFT_mass_math_stability": r"Math mass stability\s*=\s*(\S+)",
    "EFT_ghost_stability": r"\n\s*Ghost stability\s*=\s*(\S+)",
    "EFT_gradient_stability": r"Gradient stability\s*=\s*(\S+)",
    "EFT_mass_stability": r"\n\s*Mass stability\s*=\s*(\S+)",
    "EFT_mass_stability_rate": r"Mass stability rate\s*=\s*(\S+)",
    "EFT_additional_priors": r"Additional priors\s*=\s*(\S+)",
    "EFT_positivity_bounds": r"Positivity bounds\s*=\s*(\S+)",
    "EFT_minkowski_limit": r"Minkowski limit\s*=\s*(\S+)",
    "EFTflag": r"EFTflag\s*=\s*(\S+)",
    "DesignerEFTmodel": r"DesignerEFTmodel\s*=\s*(\S+)",
}


def child(tree, ini):
    sys.path.insert(0, tree)
    import camb

    cp = camb.read_ini(ini, no_validate=True)
    for k, v in VALUES.items():
        setattr(cp.EFTCAMB, k, v)
    sys.stdout.flush()
    cp.EFTCAMB.feedback(print_params=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
    ap.add_argument("--child", nargs=2)
    a = ap.parse_args()
    if a.child:
        child(*a.child)
        return
    tmp = tempfile.mkdtemp(prefix="echo_")
    ini = os.path.join(tmp, "fR.ini")
    Path(ini).write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\noutput_root = {tmp}/fR\n"
                         "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    out = subprocess.run([sys.executable, __file__, "--child", a.tree, ini], capture_output=True, text=True,
                         env=env, cwd=tmp).stdout
    expected = dict(VALUES, EFTflag=3, DesignerEFTmodel=1)
    bad, seen = {}, 0
    for k, pat in PATTERNS.items():
        m = re.search(pat, out)
        if not m:
            bad[k] = "not printed"
            continue
        seen += 1
        tok, exp = m.group(1), expected[k]
        if isinstance(exp, bool):
            ok = tok == ("T" if exp else "F")
        else:
            ok = abs(float(tok) - exp) <= 1e-5 * abs(exp)
        if not ok:
            bad[k] = [tok, exp]
    print(json.dumps({"fields_checked": seen, "mismatches": bad}, indent=1))


if __name__ == "__main__":
    main()
