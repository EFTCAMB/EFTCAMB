"""Independent check of finding eftcamb-import-fails-non-gfortran (verifier).

Instead of flipping baseconfig.gfortran by hand, build a library that genuinely uses the ifort symbol convention
(module procedures/variables exported as ``module_mp_name_``; bind(c) names verbatim, as every compiler does) by
relinking A0's release objects with ld64 aliases, and hiding the gfortran ``__module_MOD_name`` exports.
Then, in fresh subprocesses:
  orig : A0's python package + ifort-convention library   -> expect baseconfig.compiler == "ifort"; does `import camb` work?
  fixed: same, with camb/eftcamb.py always using the bind(c) name -> does `import camb` work and does a model run?

Usage: python ifort_naming_sim.py [step ...]   steps: link, import, run (default: all)
"""

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

SCR = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
A0 = SCR / "trees" / "A0"
OBJ = SCR / "build" / "A0" / "camb" / "Releaselib"
FORUTILS = SCR / "build" / "A0" / "forutils" / "Releaselib"
WORK = SCR / "work" / "verify" / "ifortsim"
LIB = WORK / "lib" / "camblib.so"
PY = sys.executable


def link() -> None:
    (WORK / "lib").mkdir(parents=True, exist_ok=True)
    syms = subprocess.run(["nm", "-gU", str(A0 / "camb" / "camblib.so")], capture_output=True, text=True).stdout
    names = [ln.split()[-1] for ln in syms.splitlines() if ln.split()]
    aliases = []
    for n in names:
        m = re.fullmatch(r"___(\w+?)_MOD_(\w+)", n)
        if m:
            aliases.append(f"{n} _{m.group(1)}_mp_{m.group(2)}_")
    (WORK / "lib" / "aliases.txt").write_text("\n".join(aliases) + "\n")
    (WORK / "lib" / "unexport.txt").write_text("___*_MOD_*\n")
    objs = sorted(str(p) for p in OBJ.glob("*.o"))
    cmd = (
        ["gfortran", "-dynamiclib", "-O3", "-fopenmp"]
        + objs
        + [f"-L{FORUTILS}", "-lforutils", "-lblas", "-llapack"]
        + [f"-Wl,-alias_list,{WORK / 'lib' / 'aliases.txt'}", f"-Wl,-unexported_symbols_list,{WORK / 'lib' / 'unexport.txt'}"]
        + ["-o", str(LIB)]
    )
    with harness.slot(1):
        r = subprocess.run(cmd, capture_output=True, text=True)
    print("link rc", r.returncode, r.stderr[-2000:])
    out = subprocess.run(["nm", "-gU", str(LIB)], capture_output=True, text=True).stdout
    exp = [ln.split()[-1] for ln in out.splitlines() if ln.split()]
    print("aliases requested:", len(aliases))
    print("exported total:", len(exp), " with _MOD_:", sum("_MOD_" in e for e in exp), " with _mp_:", sum("_mp_" in e for e in exp))
    print("init symbols exported:", [e for e in exp if "initialize_eftcamb" in e.lower()])

    for variant in ("orig", "fixed"):
        tree = WORK / variant
        if tree.exists():
            shutil.rmtree(tree)
        shutil.copytree(A0 / "camb", tree / "camb", ignore=shutil.ignore_patterns("__pycache__", "camblib.so"))
        shutil.copy2(LIB, tree / "camb" / "camblib.so")
        (tree / "AUDIT_TREE.json").write_text(json.dumps({"note": "A0 python + ifort-convention relink", "variant": variant}))
    f = WORK / "fixed" / "camb" / "eftcamb.py"
    src = f.read_text()
    old = (
        "if gfortran:\n    EFTCAMB_initialize_read_params = camblib.initialize_EFTCAMB_from_py_read_params_\n"
        "else:\n    EFTCAMB_initialize_read_params = camblib.initialize_EFTCAMB_from_py_read_params\n"
    )
    assert old in src
    f.write_text(src.replace(old, "EFTCAMB_initialize_read_params = camblib.initialize_EFTCAMB_from_py_read_params_\n"))


PROBE = r"""
import sys, traceback
sys.path.insert(0, sys.argv[1])
try:
    import camb
    from camb import baseconfig
    print("IMPORT OK compiler=%s gfortran=%s file=%s" % (baseconfig.compiler, baseconfig.gfortran, camb.__file__))
except BaseException as e:
    tb = traceback.extract_tb(e.__traceback__)
    bc = sys.modules.get("camb.baseconfig")
    print("IMPORT FAILED compiler=%s gfortran=%s" % (bc.compiler, bc.gfortran))
    print("  %s: %s" % (type(e).__name__, str(e).splitlines()[-1][:140]))
    print("  chain:", " <- ".join("%s:%d" % (f.filename.split('/camb/')[-1], f.lineno) for f in reversed(tb) if '/camb/' in f.filename))
"""


def do_import() -> None:
    for variant in ("orig", "fixed"):
        env = dict(os.environ, OMP_NUM_THREADS="1")
        env.pop("PYTHONPATH", None)
        r = subprocess.run([PY, "-c", PROBE, str(WORK / variant)], capture_output=True, text=True, env=env, cwd=WORK)
        print(f"[{variant}]", r.stdout.strip(), r.stderr.strip()[-500:])
    # control: genuine gfortran library (A0) imports fine
    r = subprocess.run([PY, "-c", PROBE, str(A0)], capture_output=True, text=True, cwd=WORK)
    print("[A0 control]", r.stdout.strip())


def run() -> None:
    ini = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
    out = WORK / "run_fixed_2_PEFT_gamma1_const_1.npz"
    st = harness.run_model(WORK / "fixed", ini, out, threads=1)
    print("fixed-tree run:", st["status"], st.get("log_tail", "")[-5:] if st["status"] != "ok" else "")
    if st["status"] == "ok":
        ref = SCR / "runs" / "A0" / "2_PEFT_gamma1_const_1.npz"
        c = harness.compare(harness.load(out), harness.load(ref))
        print("vs A0 reference: bit_identical", c["bit_identical"], "worst_cl", c["worst_cl"], "worst_pk", c["worst_pk"])
        meta = harness.load(out)["meta"]
        meta = meta if isinstance(meta, dict) else json.loads(str(meta))
        print("meta model:", meta.get("eft_model_name"), "EFTflag", meta.get("EFTflag"), "init_fallback", meta.get("eft_init_fallback"))


if __name__ == "__main__":
    steps = sys.argv[1:] or ["link", "import", "run"]
    for s in steps:
        print(f"=== {s}")
        {"link": link, "import": do_import, "run": run}[s]()
