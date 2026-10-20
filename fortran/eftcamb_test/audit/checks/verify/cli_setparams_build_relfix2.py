"""Release copy of the frozen A0 tree with the alternative fix: drop the EFTCAMB exemption and restore the upstream
CAMB check `if (.not. allocated(P%DarkEnergy)) call GlobalError('DarkEnergy not set', ...)`."""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from cli_setparams_build_debug import OLD, OUT, S  # noqa: E402

NEW2 = """    if (.not. allocated(P%DarkEnergy)) then
        call GlobalError('DarkEnergy not set', error_darkenergy)
    end if
"""
dest = OUT / "relfix2"
if dest.exists():
    shutil.rmtree(dest)
shutil.copytree(S / "trees" / "A0", dest, symlinks=True)
res = dest / "fortran" / "results.f90"
txt = res.read_text()
assert txt.count(OLD) == 1
res.write_text(txt.replace(OLD, NEW2))
so = dest / "camb" / "camblib.so"
t0 = so.stat().st_mtime
for i in range(4):
    r = harness.build_tree(S / "trees" / "A0", dest)
    log = Path(r["log"]).read_text()
    updated = so.stat().st_mtime > t0
    print(i, r["returncode"], r["elapsed"], "results.f90 compiled:", "results.f90" in log, "so updated:", updated)
    if r["returncode"] == 0 and updated:
        break
