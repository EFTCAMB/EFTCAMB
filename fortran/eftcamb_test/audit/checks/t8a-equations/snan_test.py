"""Uninitialised-variable reads in the EFT blocks of equations.f90, detected with the PS tree (code under audit +
inert audit patches, built with -O0 -finit-real=snan -ffpe-trap=invalid,zero,overflow, no -fcheck so F-001/F-002
do not interfere).

Runs 2_PEFT_Omega_const_1 (scalars + tensors, reduced lmax) three times:
  flags off                     -> expected SIGFPE on w_dark_energy_t (derivs output block, "standard DE perturbations")
  AUDIT_FIX_WDE=1               -> expected SIGFPE on grhog_t/grhor_t in derivst (EFT tensor block)
  AUDIT_FIX_WDE=1 FIX_TENSOR=1  -> expected to run (or to trap at the next uninitialised read)
and prints status and the first equations.f90 frame of each backtrace.

Usage: python snan_test.py [MODEL]
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8a-equations")
model = sys.argv[1] if len(sys.argv) > 1 else "2_PEFT_Omega_const_1"
ini = harness.PARAMS_DIR / f"{model}.ini"
ov = {"l_max_scalar": "600", "l_max_tensor": "400", "k_eta_max_tensor": "800", "get_transfer": "F",
      "do_lensing": "F", "num_redshiftwindows": "0"}
out = W / "snan"
out.mkdir(exist_ok=True)
patched = (W / "PS" / "fortran" / "equations.f90").read_text().splitlines()
res = {}
for tag, env in (("off", {}), ("wde", {"AUDIT_FIX_WDE": "1"}), ("wde_tensor", {"AUDIT_FIX_WDE": "1",
                                                                             "AUDIT_FIX_TENSOR": "1"})):
    r = harness.run_model(W / "PS", ini, out / f"{model}_{tag}.npz", overrides=ov, extra_env=env, threads=1,
                          timeout=3000)
    log = Path(r["log"]).read_text(errors="replace")
    frames = re.findall(r"(\S*equations\.f90):(\d+)", log)
    first = None
    if frames:
        ln = int(frames[0][1])
        first = {"patched_line": ln, "code": patched[ln - 1].strip()}
    res[tag] = {"status": r["status"], "first_equations_frame": first,
                "fpe": "Floating-point exception" in log or "SIGFPE" in log}
    print(tag, json.dumps(res[tag]))
(out / f"{model}.json").write_text(json.dumps(res, indent=1))
