"""Build two A0 copies: P = Kf(2) added only to the pigddot octupole term, E = Kf(2) only on E2ddot's Edot(3)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "polter"
W.mkdir(parents=True, exist_ok=True)
OLD_P = "- 3._dl/5._dl*k*octgdot - dopacity*(pig - polter)"
NEW_P = "- 3._dl/5._dl*k*EV%Kf(2)*octgdot - dopacity*(pig - polter)"
OLD_E = "opacity*(Edot(2) - polterdot) - k/3._dl*Edot(3)"
NEW_E = "opacity*(Edot(2) - polterdot) - k/3._dl*EV%Kf(2)*Edot(3)"
for name, (old, new) in {"P": (OLD_P, NEW_P), "E": (OLD_E, NEW_E)}.items():
    if name not in sys.argv[1:]:
        continue
    dest = W / name
    if not dest.exists():
        import shutil
        shutil.copytree(harness.tree_path("A0"), dest, symlinks=True)
    f = dest / "fortran" / "equations.f90"
    s = f.read_text()
    if new not in s:
        assert s.count(old) == 1, name
        f.write_text(s.replace(old, new))
    print(name, harness.build_tree("A0", dest))
