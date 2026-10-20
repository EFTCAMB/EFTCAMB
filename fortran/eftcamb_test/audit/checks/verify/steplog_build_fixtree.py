"""Adversarial verification (steplog-d4-and-integral): build a copy of A0 with ONLY the steplog d4/integral fix.

Copy A0 -> work/verify/steplog/fixtree, replace the bodies of steplogParametrized1DFourthDerivative and
steplogParametrized1DIntegral with the closed forms (independently written), build with release flags.
Pass --d4-only or --int-only to build a tree with only one of the two fixes (tree name suffix _d4 / _int).
"""

import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "steplog_phys"
REL = "fortran/eftcamb/04f_parametrizations_1D/04p018_steplog_parametrizations_1D.f90"

D4_NEW = """        term1 = (log(x) - log(self%aT))/self%Delta
        steplogParametrized1DFourthDerivative = 0.5_dl*(self%v2 - self%v1)/x**4*( &
            15._dl*term1*(3._dl - 4._dl*term1**2)/(1._dl + term1**2)**4.5_dl/self%Delta**4 &
            - 6._dl*(12._dl*term1**2 - 3._dl)/(1._dl + term1**2)**3.5_dl/self%Delta**3 &
            - 33._dl*term1/(1._dl + term1**2)**2.5_dl/self%Delta**2 &
            - 6._dl/(1._dl + term1**2)**1.5_dl/self%Delta)
"""

INT_NEW = """        term1 = x**(-1._dl - 1.5_dl*(self%v1 + self%v2))
        term2 = exp(-1.5_dl*(self%v2 - self%v1)* &
            (sqrt(self%Delta**2 + (log(x) - log(self%aT))**2) - sqrt(self%Delta**2 + log(self%aT)**2)))
"""


def main():
    d4 = "--int-only" not in sys.argv
    it = "--d4-only" not in sys.argv
    name = "fixtree" + ("" if d4 and it else ("_d4" if d4 else "_int"))
    dest = WORK / name
    WORK.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        shutil.copytree(harness.tree_path("A0"), dest, symlinks=True)
    f = dest / REL
    src = f.read_text()
    if d4:
        pat = re.compile(r"(        steplogParametrized1DFourthDerivative = \(1\._dl/2\._dl\).*?/self%Delta\n)", re.S)
        src, n = pat.subn(D4_NEW, src)
        print("d4 hunk replaced:", n)
    if it:
        pat = re.compile(r"(        term1 = x\*\*\(-1\._dl - 1\.5_dl\*\(self%v1 \+ self%v2\) - &.*?log\(self%aT\)\*\*2\)\)\)\n)",
                         re.S)
        src, n = pat.subn(INT_NEW, src)
        print("integral hunk replaced:", n)
    f.write_text(src)
    print(harness.build_tree(dest, dest))


if __name__ == "__main__":
    main()
