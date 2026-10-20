"""Repro: EFT_IC_type = 2/3 read Gamma1..3 (and the background) from a stale / uninitialised EV%eft_cache.

equations.f90 (subroutine initial, IC types 2/3) calls only
    CP%EFTCAMB%model%compute_background_EFT_functions(a_eval, CP%eft_par_cache, EV%eft_cache)
and then reads EV%eft_cache%EFTGamma1V/2V/3V.  The Gammas are set by compute_secondorder_EFT_functions, never called
here, and EV%eft_cache is not initialised at the start of a k mode (EV is OMP-private and reused), so the first mode
on each thread sees uninitialised memory and later modes see the values left by the previous mode at a ~ 1.  Pure-EFT
background functions also read eft_cache%adotoa / densities that are stale here.

The frozen release tree cannot print the values, so this script uses the audit copy
work/t5-robust/A0_icdbg, which is A0 plus (a) a print of q, a_eval, Omega0, Gamma1, Gamma2, Gamma3 at that point when
EFT_T5_DEBUG=1 and (b) a candidate fix enabled by EFT_T5_ICFIX=1 (fill the full cache at a_eval with the same call
sequence used by the RGR module).  It runs 2_PEFT_gamma3_exp_1 with EFT_IC_type = 3 on one OpenMP thread.

Expected: Gamma3 = Gamma3(a_eval = 0.01) = 0.3709 for every k mode (the value the RGR feedback prints).
Run: /Users/marco/Software/miniforge3/bin/python ic23-stale-cache-gammas.py
"""

import os
import subprocess
import sys
import tempfile

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = ROOT + "/work/t5-robust/A0_icdbg"
PARAMS = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
CHILD = "import sys; sys.path.insert(0, sys.argv[1]); import camb; cp = camb.read_ini(sys.argv[2], no_validate=True)\n" \
        "try:\n    camb.get_results(cp); print('RUN OK')\nexcept Exception as e:\n    print('CAMBError', str(e)[:160])\n"


def run(fix: str) -> list[str]:
    tmp = tempfile.mkdtemp()
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/2_PEFT_gamma3_exp_1.ini)\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                "EFTCAMB_write_background = F\nEFT_IC_type = 3\n")
    env = dict(os.environ, OMP_NUM_THREADS="1", EFT_T5_DEBUG="1", EFT_T5_ICFIX=fix)
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini], capture_output=True, text=True, cwd=tmp, env=env)
    lines = (p.stdout + p.stderr).splitlines()
    return [x for x in lines if x.startswith("T5IC") or x.startswith("RUN OK") or x.startswith("CAMBError")]


for fix in ("0", "1"):
    print(f"--- EFT_T5_ICFIX={fix} ({'release code path' if fix == '0' else 'candidate fix'})")
    for line in run(fix)[:6]:
        parts = line.split()
        if line.startswith("T5IC"):
            q, a_eval, om, g1, g2, g3 = parts[-6:]
            print(f"  k-mode q={float(q):9.4f} a_eval={float(a_eval):.3g} Omega0={om} Gamma1={g1} Gamma2={g2} "
                  f"Gamma3={g3}   (expected Gamma1=0 Gamma2=0 Gamma3=0.3709)")
        else:
            print("  " + line)
