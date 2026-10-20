"""Build A0 copy with only -finit-real=snan -finit-derived added (no source edits)."""
import sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
W = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/uninit"
COMMON = "-MMD -cpp -ffree-line-length-none -fallow-argument-mismatch -fmax-errors=4 -fopenmp -lblas -llapack"
print(harness.build_tree(sys.argv[1] if len(sys.argv)>1 else "A0", f"{W}/snanfix", fflags=f"-fPIC -O3 -w -finit-real=snan -finit-derived {COMMON}"))
