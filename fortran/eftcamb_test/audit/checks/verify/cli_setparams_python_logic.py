"""Second half of the finding: the 'DarkEnergy not set' test reads this%CP (the *previous* parameters of the CAMBdata)
instead of the new P. Exercise it from Python (where CAMBdata.__init__ allocates this%CP%EFTCAMB, so no -O0 crash)
on the frozen release A0 and on a release copy with only the proposed fix (relfix).

Cases (each in its own subprocess; a segfault only kills that case):
  T0 fresh CAMBdata, 4_Des_F_R_LCDM_1 (EFT on), DarkEnergy allocated           -> sanity, sigma8
  T1 fresh CAMBdata, 4_Des_F_R_LCDM_1 (EFT on), pars.DarkEnergy = None
  T2 fresh CAMBdata, 1_EFT_GR (EFTflag=0),      pars.DarkEnergy = None
  T3 CAMBdata reused after a T0 run, then 1_EFT_GR (EFTflag=0) with pars.DarkEnergy = None

relfix2 = release copy with the upstream check restored (no EFTCAMB exemption).

Usage: python cli_setparams_python_logic.py [A0] [relfix] [relfix2]
"""

import json
import os
import signal
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

PAR = Path("/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters")
S = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
    "/scratchpad/audit"
)
TREES = {"A0": S / "trees" / "A0", "relfix": S / "work" / "verify" / "cli_setparams" / "relfix",
         "relfix2": S / "work" / "verify" / "cli_setparams" / "relfix2"}
if len(sys.argv) > 1:
    TREES = {k: v for k, v in TREES.items() if k in sys.argv[1:]}
OUTDIR = S / "work" / "verify" / "cli_setparams" / "pylogic"

CHILD = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
import camb
case = sys.argv[2]
des, gr = sys.argv[3], sys.argv[4]
def load(ini):
    p = camb.read_ini(ini)
    return p
res = {"case": case}
try:
    if case == "T0":
        p = load(des)
        r = camb.get_results(p)
        res["sigma8"] = float(r.get_sigma8_0())
        res["EFTflag"] = int(p.EFTCAMB.EFTflag)
    elif case in ("T1", "T2"):
        p = load(des if case == "T1" else gr)
        p.DarkEnergy = None
        res["DarkEnergy_after_None"] = repr(p.DarkEnergy)
        res["EFTflag"] = int(p.EFTCAMB.EFTflag)
        r = camb.get_results(p)
        res["sigma8"] = float(r.get_sigma8_0())
    elif case == "T3":
        p1 = load(des)
        r = camb.get_results(p1)
        res["first_sigma8"] = float(r.get_sigma8_0())
        p2 = load(gr)
        p2.DarkEnergy = None
        res["EFTflag_second"] = int(p2.EFTCAMB.EFTflag)
        r.calc_power_spectra(p2)
        res["sigma8"] = float(r.get_sigma8_0())
    res["outcome"] = "ok"
except Exception as e:
    res["outcome"] = f"exception {type(e).__name__}: {e}"
print("RESULT " + json.dumps(res), flush=True)
"""


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    child = OUTDIR / "child.py"
    child.write_text(CHILD)
    for tname, tree in TREES.items():
        for case in ("T0", "T1", "T2", "T3"):
            with harness.slot(1):
                p = subprocess.Popen(
                    [sys.executable, "-I", str(child), str(tree), case, str(PAR / "4_Des_F_R_LCDM_1.ini"),
                     str(PAR / "1_EFT_GR.ini")], cwd=OUTDIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, start_new_session=True, env=dict(os.environ, OMP_NUM_THREADS="1"),
                    preexec_fn=lambda: os.nice(10))
                try:
                    out, _ = p.communicate(timeout=400)
                    rc = p.returncode
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                    out, _ = p.communicate()
                    rc = None
            (OUTDIR / f"{tname}_{case}.log").write_text(out)
            line = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
            if line:
                summary = json.loads(line[0][7:])
            elif rc is None:
                summary = {"outcome": "timeout"}
            elif rc < 0:
                summary = {"outcome": f"killed by {signal.Signals(-rc).name}"}
            else:
                summary = {"outcome": f"rc={rc}", "tail": out.splitlines()[-3:]}
            print(f"{tname:7s} {case}: {json.dumps(summary)}")


if __name__ == "__main__":
    main()
