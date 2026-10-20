"""Verify ic23-stale-cache-gammas: (1) rerun instrumented tree with correct paths; (2) unmodified A0 thread dependence."""
import os, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness
import re
fl = lambda v: float(re.sub(r'(\d)([-+]\d{3})$', r'\1E\2', v))
import numpy as np

S = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
OUT = Path(S) / "work/verify/ic23"
PARAMS = harness.REPO / "fortran/eftcamb_test/parameters"
mode = sys.argv[1]
if mode == "instr":
    TREE = S + "/work/t5-robust/A0_icdbg"
    CHILD = ("import sys; sys.path.insert(0, sys.argv[1]); import camb; cp = camb.read_ini(sys.argv[2], no_validate=True)\n"
             "camb.get_results(cp); print('RUN OK')\n")
    for model in sys.argv[2:]:
        for fix in ("0", "1"):
            tmp = tempfile.mkdtemp(dir=OUT)
            ini = os.path.join(tmp, "w.ini")
            Path(ini).write_text(f"DEFAULT({PARAMS}/{model}.ini)\noutput_root = {tmp}/r\nfeedback_level = 1\n"
                                 "EFTCAMB_write_background = F\nEFT_IC_type = 3\n")
            env = dict(os.environ, OMP_NUM_THREADS="1", EFT_T5_DEBUG="1", EFT_T5_ICFIX=fix)
            env.pop("PYTHONPATH", None)
            with harness.slot(1):
                p = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, TREE, ini], capture_output=True,
                                   text=True, cwd=tmp, env=env, timeout=1500)
            lines = [x for x in (p.stdout + p.stderr).splitlines() if x.startswith("T5IC")]
            print(model, "fix", fix, "rc", p.returncode, "nmodes", len(lines)); print("   tail:", (p.stdout+p.stderr).strip().splitlines()[-3:])
            g3 = np.array([[fl(v) for v in l.split()[-6:]] for l in lines])
            for r in g3[:4]:
                print("   q=%.4g a=%.3g Om=%.5g G1=%.4g G2=%.4g G3=%.5g" % tuple(r))
            if len(g3):
                print("   unique G3 values:", np.unique(g3[:, 5])[:8], "count", len(np.unique(g3[:, 5])))
else:
    model = sys.argv[2]
    ini = PARAMS / f"{model}.ini"
    res = {}
    for th in (1, 2, 4):
        o = OUT / f"A0_{model}_ic3_t{th}.npz"
        st = harness.run_model("A0", ini, o, overrides={"EFT_IC_type": "3"}, threads=th)
        print(th, st["status"])
        res[th] = o
    for th in (2, 4):
        c = harness.compare(res[1], res[th])
        print(f"t1 vs t{th}: worst_cl", c["worst_cl"], "bit_identical", c["bit_identical"])
    o = OUT / f"A0_{model}_ic3_t1b.npz"
    harness.run_model("A0", ini, o, overrides={"EFT_IC_type": "3"}, threads=1)
    c = harness.compare(res[1], o); print("t1 vs t1 rerun: worst_cl", c["worst_cl"], "bit", c["bit_identical"])
