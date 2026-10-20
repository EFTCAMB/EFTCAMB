"""Attribution check (verify) for initialize-parameters-bad-value-kills-python: run the python repro case
EFTflag=1.0 under lldb (batch mode), break on gfortran's error-stop entry points and print the call stack, to confirm
the process is stopped by forutils TNameValueList_FailStop reached from the EFTCAMB python bridge Ini.

Uses the AD tree (debug symbols, -O0) so frames carry file:line; F-001/F-002 code paths are not reached here
(no character getters, no thermo init).

Run: /Users/marco/Software/miniforge3/bin/python initparams_badvalue_lldb.py [--tree A0|AD]
"""

import argparse
import os
import subprocess
import sys

TREES = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees"
WORK = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/work/verify/initparams_badvalue")

CHILD = """import sys
sys.path.insert(0, {tree!r})
import camb
p = camb.CAMBparams()
p.set_cosmology(H0=67.5)
d = {{'EFTflag': 1.0, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1, 'EFTOmega0': 0.01, 'EFTwDE': 0}}
p.EFTCAMB.initialize_parameters(p, d, print_header=False)
print('CHILD survived', flush=True)
"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="AD")
    a = ap.parse_args()
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import harness

    os.makedirs(WORK, exist_ok=True)
    script = os.path.join(WORK, f"child_{a.tree}.py")
    with open(script, "w") as f:
        f.write(CHILD.format(tree=f"{TREES}/{a.tree}"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    cmd = ["lldb", "--batch",
           "-o", "settings set target.process.stop-on-exec false",
           "-o", "breakpoint set -r '^_gfortran_error_stop'",
           "-o", "breakpoint set -n exit",
           "-o", "run",
           "-o", "bt 12",
           "-k", "bt 12", "-k", "kill",
           "--", sys.executable, script]
    with harness.slot(1):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=500)
    out = r.stdout + r.stderr
    with open(os.path.join(WORK, f"lldb_{a.tree}.log"), "w") as f:
        f.write(out)
    keep = [ln for ln in out.splitlines()
            if "frame #" in ln or "Error for key" in ln or "stop reason" in ln or "CHILD" in ln or "Breakpoint" in ln]
    print("\n".join(keep[:60]))
