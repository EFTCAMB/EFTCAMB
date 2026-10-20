"""Build and run interp_xfinal_unit.f90 against a frozen tree's objects (AD = -fcheck=all; A0 = release)."""
import glob, os, subprocess, sys
AUD = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
tree = sys.argv[1] if len(sys.argv) > 1 else "AD"
L = f"{AUD}/build/{tree}/camb/Releaselib"
F = f"{AUD}/build/{tree}/forutils/Releaselib"
W = f"{AUD}/work/verify"
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "interp_xfinal_unit.f90")
objs = [o for o in sorted(glob.glob(f"{L}/0[12]*.o")) if "python" not in o] + [f"{L}/{n}.o" for n in ("constants", "config", "MathUtils")]
exe = f"{W}/interp_unit_{tree}"
flags = ["-g", "-O0", "-fcheck=all"] if tree.startswith("AD") else ["-O3"]
r = subprocess.run(["gfortran", *flags, "-cpp", "-ffree-line-length-none", "-fopenmp", f"-I{L}", "-J", W, src, *objs, f"{F}/libforutils.a", "-o", exe, "-lblas", "-llapack"], capture_output=True, text=True)
errs = [l for l in r.stderr.splitlines() if "warning" not in l]
if r.returncode:
    print("\n".join(errs[-15:])); sys.exit(1)
r = subprocess.run([exe], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
print(r.stdout + r.stderr, "rc=", r.returncode)
