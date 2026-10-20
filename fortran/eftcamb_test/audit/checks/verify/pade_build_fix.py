"""Build a copy of A0 with only the one-line Pade fourth-derivative fix (finding pade-fourth-derivative)."""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify")
dest = W / "pade_fix"
src = harness.tree_path("A0")
if not dest.exists():
    shutil.copytree(src, dest, symlinks=True)
f = dest / "fortran/eftcamb/04f_parametrizations_1D/04p010_pade_series_1D.f90"
s = f.read_text()
old = "        4._dl*(-ddddown/down**2 - 6._dl*ddown**3/down**4 + 6._dl*ddown*dddown/down**3) + &"
new = "        4._dl*dup*(-ddddown/down**2 - 6._dl*ddown**3/down**4 + 6._dl*ddown*dddown/down**3) + &"
if old in s:
    assert s.count(old) == 1
    f.write_text(s.replace(old, new))
assert new in f.read_text()
print(harness.build_tree(src, dest))
