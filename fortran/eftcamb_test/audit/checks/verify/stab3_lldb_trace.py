"""Independent check of finding stability-third-pass-stale-y: trace, with lldb breakpoints, the scale factor passed to
every EFTTestStability call made by EFTCAMB_Stability_Check in an unmodified build (no DebugEFTCAMB dump involved).

The breakpoint callback reads the double pointed to by x0 (first by-reference argument `a` of EFTTestStability, arm64
ABI) and the astart/aend arguments (x3, x4) of EFTCAMB_Stability_Check.

Usage: python stab3_lldb_trace.py [TREE] [INI] [OUTDIR]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
             "/scratchpad/audit")
REPO = Path(__file__).resolve().parents[5]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else AUDIT / "trees/A0"
INI = Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / "fortran/eftcamb_test/parameters/2_PEFT_Omega_const_1.ini"
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else AUDIT / "work/verify/stab3"
OUT.mkdir(parents=True, exist_ok=True)

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
camb.get_background(pars)
print("CHILD_DONE", flush=True)
"""

CALLBACKS = r'''
import lldb, struct
LOG = open(__LOGPATH__, "w")

def _dbl(frame, reg):
    proc = frame.GetThread().GetProcess()
    addr = frame.FindRegister(reg).GetValueAsUnsigned()
    err = lldb.SBError()
    raw = proc.ReadMemory(addr, 8, err)
    return struct.unpack("<d", raw)[0] if err.Success() else float("nan")

def on_test(frame, bp_loc, extra_args, internal_dict):
    LOG.write("T %.17e\n" % _dbl(frame, "x0"))
    return False

def on_check(frame, bp_loc, extra_args, internal_dict):
    LOG.write("C %.17e %.17e\n" % (_dbl(frame, "x3"), _dbl(frame, "x4")))
    LOG.flush()
    return False

def on_exit_flush(*a):
    LOG.flush()
'''


def main():
    tmp = Path(tempfile.mkdtemp(prefix="stab3lldb_", dir=OUT))
    log = tmp / "trace.txt"
    cb = tmp / "cb.py"
    cb.write_text(CALLBACKS.replace("__LOGPATH__", repr(str(log))))
    child = tmp / "child.py"
    child.write_text(CHILD)
    w = tmp / "w.ini"
    w.write_text(f"DEFAULT({INI})\noutput_root = {tmp / 'run'}\nfeedback_level = 0\nEFTCAMB_write_background = F\n")
    cmds = [
        f"command script import {cb}",
        "settings set target.process.stop-on-exec false",
        "breakpoint set -n __eftcamb_stability_MOD_eftteststability",
        "breakpoint command add -F cb.on_test 1",
        "breakpoint set -n __eftcamb_stability_MOD_eftcamb_stability_check",
        "breakpoint command add -F cb.on_check 2",
        "run",
        "script cb.on_exit_flush()",
        "breakpoint list",
        "quit",
    ]
    argv = ["lldb", "-b"]
    for c in cmds:
        argv += ["-o", c]
    argv += ["--", sys.executable, str(child), str(TREE), str(w)]
    with harness.slot(1):
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        r = subprocess.run(argv, capture_output=True, text=True, cwd=tmp, env=env, timeout=540,
                           preexec_fn=lambda: os.nice(10))
    (tmp / "lldb_stdout.txt").write_text(r.stdout + "\n---stderr---\n" + r.stderr)
    lines = log.read_text().splitlines() if log.exists() else []
    checks = [tuple(map(float, l.split()[1:])) for l in lines if l.startswith("C ")]
    tests = np.array([float(l.split()[1]) for l in lines if l.startswith("T ")])
    print(f"tree: {TREE}\nini: {INI.name}\nlldb rc={r.returncode}, CHILD_DONE={'CHILD_DONE' in r.stdout}")
    print(f"EFTCAMB_Stability_Check calls: {len(checks)}  (astart, aend) = {checks}")
    print(f"EFTTestStability calls traced: {tests.size}")
    res = {"tree": str(TREE), "ini": str(INI), "n_check_calls": len(checks), "check_args": checks,
           "n_test_calls": int(tests.size), "trace": str(log)}
    if tests.size >= 3000 and checks:
        astart, aend = checks[0]
        passes = {"pass1": tests[:1000], "pass2": tests[1000:2000], "pass3": tests[2000:3000]}
        for name, p in passes.items():
            print(f"{name}: unique a = {np.unique(p).size:4d}, min = {p.min():.10e}, max = {p.max():.10e}, "
                  f"min(aend - a) = {np.min(aend - p):.3e}")
            res[name] = {"unique": int(np.unique(p).size), "min": float(p.min()), "max": float(p.max()),
                         "min_aend_minus_a": float(np.min(aend - p))}
        i = np.arange(1000)
        y = -10.0 + i * 10.0 / 999
        res["pass1_matches_formula"] = bool(np.allclose(passes["pass1"], astart + i * (aend - astart) / 999,
                                                        rtol=1e-14, atol=0))
        res["pass2_matches_formula"] = bool(np.allclose(passes["pass2"], astart + (aend - astart) * 10**y,
                                                        rtol=1e-12, atol=0))
        res["pass3_all_equal_astart"] = bool(np.all(passes["pass3"] == astart))
        intended = aend + (astart - aend) * 10**y
        res["pass3_intended_min_aend_minus_a"] = float(np.min(aend - intended))
        print("pass1 matches linear formula:", res["pass1_matches_formula"])
        print("pass2 matches log-near-astart formula:", res["pass2_matches_formula"])
        print("pass3: every a == astart:", res["pass3_all_equal_astart"])
        print(f"pass3 intended: min(aend - a) = {res['pass3_intended_min_aend_minus_a']:.3e}")
    else:
        print(r.stdout[-3000:], r.stderr[-2000:])
    (OUT / f"lldb_trace_{TREE.name}_{INI.stem}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
