"""Repro: the third stability sampling pass ("log sampling close to aend") never moves: it re-tests a = astart 1000 times.

Code: fortran/eftcamb/09_EFTCAMB_stability.f90:122-136
        do ind = 1, indMax
            Atest = aend + (astart - aend)*10._dl**y       <- y is not updated in this loop
The second pass leaves y = LogSamplingScale + (indMax-1)*(0-LogSamplingScale)/(indMax-1) = 0, so every Atest of the
third pass equals aend + (astart - aend) = astart.  The intended dense sampling next to aend (a -> 1) is lost.

Evidence from the unmodified debug tree ADD (DebugEFTCAMB = .true.): EFTStabilityComputation dumps every tested a to
<output_root>_stability_cache_FRW.dat (first column).  The run is stopped by the known F-002 -fcheck abort in thermo
init right after the stability check, which does not matter here (the stability files are closed before).

Usage: /Users/marco/Software/miniforge3/bin/python stability-third-pass-stale-y.py [ADD_TREE]
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/ADD"
)
INI = REPO / "fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
camb.get_results(pars)
"""


def main():
    tmp = Path(tempfile.mkdtemp(prefix="stab3_"))
    root = tmp / "run"
    w = tmp / "w.ini"
    w.write_text(f"DEFAULT({INI})\noutput_root = {root}\nfeedback_level = 1\nEFTCAMB_write_background = F\n")
    r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w)], capture_output=True, text=True, cwd=tmp,
                       env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=1200)
    f = Path(str(root) + "_stability_cache_FRW.dat")
    if not f.exists():
        print("no stability dump found; rc =", r.returncode, r.stdout[-1500:], r.stderr[-1500:])
        sys.exit(1)
    a = np.loadtxt(f, usecols=0)
    print(f"tree: {TREE}  (child rc={r.returncode}; non-zero is the known F-002 abort after the stability check)")
    print(f"stability points dumped: {a.size} (3 passes x 1000 expected)")
    if a.size < 3000:
        sys.exit(1)
    p1, p2, p3 = a[:1000], a[1000:2000], a[2000:3000]
    astart, aend = p1[0], p1[-1]
    print(f"astart = {astart:.6e}, aend = {aend:.6e}")
    print(f"pass 1 (linear):            unique a = {np.unique(p1).size:4d}, range [{p1.min():.6e}, {p1.max():.6e}]")
    print(f"pass 2 (log near astart):   unique a = {np.unique(p2).size:4d}, range [{p2.min():.6e}, {p2.max():.6e}]")
    print(f"pass 3 (log near aend):     unique a = {np.unique(p3).size:4d}, range [{p3.min():.6e}, {p3.max():.6e}]")
    y = -10.0 + np.arange(1000) * 10.0 / 999
    expect = aend + (astart - aend) * 10.0**y
    print(f"expected pass 3: unique a = {np.unique(expect).size}, range [{expect.min():.6e}, {expect.max():.6e}],"
          f" closest to aend: aend - a = {aend - expect.max():.3e}")
    print(f"observed pass 3: closest to aend: aend - a = {aend - p3.max():.3e}")


if __name__ == "__main__":
    main()
