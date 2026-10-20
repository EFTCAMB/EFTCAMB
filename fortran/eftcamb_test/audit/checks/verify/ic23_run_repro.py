"""Run the ic23-stale-cache-gammas repro under a harness CPU slot, fixing the moved repo path."""
import pathlib, subprocess, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import harness
src = pathlib.Path(__file__).resolve().parents[2] / "repro" / "ic23-stale-cache-gammas.py"
code = src.read_text().replace("/Users/marco/Desktop/IW/", "/Users/marco/IW/")
with harness.slot(1):
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600)
print(p.stdout, p.stderr[-2000:])
