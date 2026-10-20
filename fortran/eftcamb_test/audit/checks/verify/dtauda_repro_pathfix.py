"""Run the t3-stability-dtauda-without-matter repro unchanged except for the ini DEFAULT path (repo moved)."""

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
src = REPO / "fortran/eftcamb_test/audit/inis/007p5-OL-gamma.ini"
OUTD = Path(sys.argv[1])
ini = OUTD / "007p5-OL-gamma_fixed.ini"
ini.write_text(src.read_text().replace("/Users/marco/Desktop/IW/", "/Users/marco/IW/"))
spec = importlib.util.spec_from_file_location(
    "rep", REPO / "fortran/eftcamb_test/audit/repro/t3-stability-dtauda-without-matter.py")
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)
rep.INI = ini
rep.OUT = OUTD
sys.argv = [sys.argv[0]] + sys.argv[2:]
rep.main()
