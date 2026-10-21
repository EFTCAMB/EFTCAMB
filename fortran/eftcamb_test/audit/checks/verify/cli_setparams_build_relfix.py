"""Release (-O3, setup.py make) copy of the frozen A0 tree with only the proposed results.f90:365 fix."""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from cli_setparams_build_debug import NEW, OLD, OUT, S  # noqa: E402

dest = OUT / "relfix"
if dest.exists():
    shutil.rmtree(dest)
shutil.copytree(S / "trees" / "A0", dest, symlinks=True)
res = dest / "fortran" / "results.f90"
txt = res.read_text()
assert txt.count(OLD) == 1
res.write_text(txt.replace(OLD, NEW))
print(harness.build_tree(S / "trees" / "A0", dest))
