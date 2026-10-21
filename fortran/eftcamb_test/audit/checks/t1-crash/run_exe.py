"""Run the command-line camb executable of a (runtime-checked) tree on EFTCAMB test inis, from scratch cwds.

The Fortran main program enables -ffpe-trap and the gfortran backtrace handler, unlike the Python entry point.
Each model runs in OUTDIR/<model>/ with a wrapper ini (DEFAULT(model ini), output_root inside the run dir,
number_of_threads = 1). Output .dat files are scanned for NaN/Inf.

Usage: python run_exe.py EXE OUTDIR --models a,b,c [--jobs 4] [--timeout 7200] [--set key=value ...]
       python run_exe.py EXE OUTDIR --inis /abs/x.ini,/abs/y.ini
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402


def scan_dat(d: Path) -> dict:
    bad = {}
    for f in sorted(d.glob("*.dat")):
        try:
            txt = f.read_text(errors="replace")
        except OSError:
            continue
        n = len(re.findall(r"(?i)\bnan\b|\binfinity\b|\binf\b", txt))
        if n:
            bad[f.name] = n
    return bad


def run_one(exe: Path, ini: Path, outdir: Path, over: dict, timeout: float) -> dict:
    d = outdir / ini.stem
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    tree = exe.parents[1]
    tpl = tree / "fortran" / "HighLExtrapTemplate_lenspotentialCls.dat"
    if tpl.exists():
        shutil.copy(tpl, d / tpl.name)
    lines = [f"DEFAULT({ini.resolve()})", f"output_root = {d / ini.stem}", "number_of_threads = 1"]
    lines += [f"{k} = {v}" for k, v in over.items()]
    w = d / "run.ini"
    w.write_text("\n".join(lines) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    t0 = time.time()
    try:
        with harness.slot(1), open(d / "run.log", "w") as f:
            t0 = time.time()
            p = subprocess.run([str(exe), str(w)], cwd=d, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=timeout,
                               preexec_fn=harness._nice)
        rc = p.returncode
        status = "ok" if rc == 0 else (f"signal:{signal.Signals(-rc).name}" if rc < 0 else f"error:rc={rc}")
    except subprocess.TimeoutExpired:
        rc, status = None, "timeout"
    log = (d / "run.log").read_text(errors="replace")
    info = {"model": ini.stem, "status": status, "rc": rc, "elapsed": round(time.time() - t0, 1),
            "nonfinite_dat": scan_dat(d), "dir": str(d)}
    m = re.search(r"Program received signal.*|Fortran runtime error.*|ERROR STOP.*|EFTCAMB error.*", log)
    if m or status != "ok":
        info["first_error"] = m.group(0) if m else ""
        info["log_tail"] = log.splitlines()[-25:]
    (d / "status.json").write_text(json.dumps(info, indent=2))
    return info


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("exe")
    p.add_argument("outdir")
    p.add_argument("--models", default="")
    p.add_argument("--inis", default="")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--timeout", type=float, default=7200)
    p.add_argument("--set", nargs="*", default=[])
    a = p.parse_args()
    if a.inis:
        inis = [Path(x) for x in a.inis.split(",")]
    else:
        want = a.models.split(",")
        inis = [i for i in harness.model_inis() if i.stem in want]
    over = dict(s.split("=", 1) for s in a.set)
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda i: run_one(Path(a.exe).resolve(), i, Path(a.outdir).resolve(), over, a.timeout), inis))
    for r in res:
        print(f"{r['model']:34s} {r['status']:16s} {r['elapsed']:8.1f}s nonfinite={r['nonfinite_dat']} "
              f"{r.get('first_error', '')}")


if __name__ == "__main__":
    main()
