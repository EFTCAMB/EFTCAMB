"""Repro: copy.deepcopy / pickle of a CAMBparams carrying an EFTCAMB model drops the model and segfaults.

CAMBparams pickling (documented in camb/model.py: "save and restore a CAMBparams instance using the repr and eval
functions, or pickle it") goes through CAMB_Structure.__getstate__/__setstate__, which only carry the ctypes _fields_.
The EFTCAMB python class maps only the flags (EFTflag, ...), not the allocatable `model` component, so the restored
object has EFTflag /= 0 and an unallocated model.  get_results / get_background then dereference
State%CP%EFTCAMB%model (results.f90:584, model%initialize_background) and the process dies with SIGSEGV.
num_params(), param_values(), ... crash the same way (11_EFTCAMB_python.f90: EFTCAMB%model%parameter_number).

Expected: deepcopy/pickle either reproduce the model (H(z=1) identical to the original) or raise a Python error.
Run: /Users/marco/Software/miniforge3/bin/python pickle-deepcopy-eftcamb-segfault.py [tree]
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
from pathlib import Path

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)


def child(tree: str, mode: str) -> None:
    sys.path.insert(0, tree)
    import copy
    import pickle

    import camb

    cp = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, EFTflag=1, PureEFTmodel=1, PureEFTmodelOmega=1,
                         EFTwDE=0, EFTOmega0=1.0, feedback_level=0)
    if mode == "original":
        c2 = cp
    elif mode == "copy":
        c2 = cp.copy()
    elif mode == "deepcopy":
        c2 = copy.deepcopy(cp)
    elif mode == "pickle":
        c2 = pickle.loads(pickle.dumps(cp))
    print(f"  {mode}: EFTflag={c2.EFTCAMB.EFTflag} model_name={c2.EFTCAMB.model_name()!r}", flush=True)
    r = camb.get_background(c2)
    print(f"  {mode}: H(z=1) = {r.hubble_parameter(1.0):.12f}", flush=True)


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else "A0"
    tree = arg if os.path.isabs(arg) else str(ROOT / "trees" / arg)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    print(f"tree: {tree}")
    bad = 0
    for mode in ("original", "copy", "deepcopy", "pickle"):
        p = subprocess.run([sys.executable, __file__, "_child", tree, mode], env=env, capture_output=True, text=True)
        out = p.stdout.strip()
        print(out)
        if p.returncode < 0:
            status = f"killed by {signal.Signals(-p.returncode).name}"
        elif p.returncode > 0:
            status = f"exit code {p.returncode} (139 = SIGSEGV via shell)\n" + p.stderr.strip()[-400:]
        else:
            status = "ok"
        if status != "ok":
            bad += 1
        print(f"  {mode}: process status = {status}")
    print("expected: all four ok with the same H(z=1); observed:", "OK" if bad == 0 else f"{bad} crashed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_child":
        child(sys.argv[2], sys.argv[3])
    else:
        main()
