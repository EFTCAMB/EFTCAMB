"""Build a copy of A0 in which only the GR 00-Einstein-equation substitution in the 21cm line_distortions source is
replaced by the kinematic form when EFTCAMB is active (diagnostic, not a proposed full fix):

    code:     s(1) + [dgrho term of s(2)] = k/(3 adotoa) (sigma - z_GR) wing exptau, z_GR = (etak + dgrho/(2k))/adotoa
    patched:  add k/(3 adotoa) (z_GR - z) wing exptau  ->  k/(3 adotoa) (sigma - z) wing exptau  (z = true EFT z)

The difference between the 21cm C_l of A0 and of this tree measures (a lower bound of) the error caused by the GR
substitutions for an EFT model.
"""

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

SCR = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
SRC = SCR / "trees" / "A0"
DEST = SCR / "work" / "verify" / "21cm" / "A0_zfix"

OLD = """                    t(0) = s(1) + s(2)

                    sources(3 + w_ix) = sources(3 + w_ix) + t(0)
                end if

                if (CP%SourceTerms%line_extra) then"""
NEW = """                    t(0) = s(1) + s(2)
                    ! AUDIT DIAGNOSTIC: replace the GR 00-equation z by the true (EFT) z in the leading term
                    if (CP%EFTCAMB%EFTflag /= 0 .and. EV%EFTCAMBactive) then
                        t(0) = t(0) + k/(3._dl*adotoa)*w%wing(j)*exptau*((etak + 0.5_dl*dgrho/k)/adotoa - z)
                    end if

                    sources(3 + w_ix) = sources(3 + w_ix) + t(0)
                end if

                if (CP%SourceTerms%line_extra) then"""

if __name__ == "__main__":
    if not DEST.exists():
        shutil.copytree(SRC, DEST, symlinks=True)
    f = DEST / "fortran" / "equations.f90"
    s = f.read_text()
    if "AUDIT DIAGNOSTIC" not in s:
        assert s.count(OLD) == 1, s.count(OLD)
        f.write_text(s.replace(OLD, NEW))
    print(harness.build_tree(SRC, DEST))
