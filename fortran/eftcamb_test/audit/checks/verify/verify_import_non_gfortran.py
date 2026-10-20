"""Independent check: first-time import of camb.model with baseconfig in a simulated non-gfortran state.
Name mangler for the fake compiler maps __mod_MOD_name -> gfortran names (so every module-procedure lookup
succeeds, as it would on a real ifort lib); only bind(c) names are looked up verbatim."""
import os, sys, types, traceback
os.environ["OMP_NUM_THREADS"] = "1"
TREE = sys.argv[1]
pkg = types.ModuleType("camb"); pkg.__path__ = [os.path.join(TREE, "camb")]; sys.modules["camb"] = pkg
import camb.baseconfig as bc
print("detected:", bc.compiler)
bc._NAME_MANGLERS["fake"] = bc._gfortran_mangle
bc.compiler = "fake"; bc.gfortran = False
try:
    import camb.model  # noqa
    print("import camb.model OK")
except Exception as e:
    tb = traceback.extract_tb(e.__traceback__)
    print("FAIL", type(e).__name__, str(e)[:120]); print("at", [(os.path.basename(f.filename), f.lineno) for f in tb][-2:])
