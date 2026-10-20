"""Repro (t9-wrapper): `import camb` fails for every library not built with gfortran (ifort, flang), because
camb/eftcamb.py looks up the bind(c) routine under a different name for non-gfortran builds.

fortran/eftcamb/11_EFTCAMB_python.f90:66-67 fixes the C name for all compilers:
    subroutine initialize_EFTCAMB_from_py_read_params(...) bind(c, name='initialize_EFTCAMB_from_py_read_params_')
camb/eftcamb.py:23-27:
    if gfortran:  EFTCAMB_initialize_read_params = camblib.initialize_EFTCAMB_from_py_read_params_
    else:         EFTCAMB_initialize_read_params = camblib.initialize_EFTCAMB_from_py_read_params   <- no such symbol
CAMB v2 detects ifort/flang libraries (camb/baseconfig.py _NAME_MANGLERS, sets baseconfig.gfortran = False) and
remaps only "__module_MOD_name" lookups; a bind(c) name is exported verbatim by every compiler.  camb/model.py:29
imports camb.eftcamb, so the AttributeError aborts `import camb` itself.

No ifort here: the repro (1) lists the exported symbol of the gfortran library, (2) re-executes camb/eftcamb.py's
module code with baseconfig.gfortran = False, i.e. exactly the branch an ifort/flang build takes.
Expected: the lookup succeeds (same symbol name for every compiler).  Observed: AttributeError.

Run: /Users/marco/Software/miniforge3/bin/python eftcamb-import-fails-non-gfortran.py [TREE]
"""

import ctypes
import importlib
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0")
sys.path.insert(0, TREE)
import camb  # noqa: E402
from camb import baseconfig  # noqa: E402

lib = ctypes.CDLL(os.path.join(TREE, "camb", "camblib.so"))
for name in ("initialize_EFTCAMB_from_py_read_params_", "initialize_EFTCAMB_from_py_read_params"):
    print(f"symbol {name!r:45s} exported: {hasattr(lib, name)}")

print("detected compiler:", baseconfig.compiler)
baseconfig.gfortran = False  # what baseconfig sets for an ifort / flang library
try:
    importlib.reload(camb.eftcamb)
    print("non-gfortran branch of camb/eftcamb.py: import OK (expected)")
except AttributeError as e:
    print(f"non-gfortran branch of camb/eftcamb.py: AttributeError: {str(e).splitlines()[-1][:110]}")
    print("expected: import OK; observed: AttributeError -> `import camb` fails for ifort/flang builds")
