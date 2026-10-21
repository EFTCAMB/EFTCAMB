"""Patch a copy of the B1 tree (pre-rebase EFTCAMB) so that ThermalNuBackground%pidot/pidotdot/pidotdotdot use the
A0 (code under audit) implementation: analytic log(a m) derivatives of the CAMB v2 pressure fit plus the exact
chain rule.  Used to prove that A0-B1 differences in massive-neutrino EFT models come from the B1 derivative
functions (table derivatives with up to 40%/200% errors in d2P/d3P, and a dimensionally wrong pidotdotdot).

Usage: python patch_b1_nufix.py B1_COPY_DIR A0_TREE_DIR
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def main() -> None:
    b1, a0 = Path(sys.argv[1]), Path(sys.argv[2])
    fb = b1 / "fortran" / "massive_neutrinos.f90"
    sb = fb.read_text()
    sa = (a0 / "fortran" / "massive_neutrinos.f90").read_text()
    # A0 fit constants needed by the derivative routine
    consts = []
    for name in ("nu_rhop_am_min", "nu_rhop_am_max", "nu_rhop_fit_c", "nu_rhop_fit_inv_zmax", "nu_fit_p_denom_scale"):
        m = re.search(rf"^\s*real\(dl\), parameter :: {name} = .*?(?:\n(?:\s+.*&\n)*\s+.*)?$", sa, re.M)
        consts.append(m.group(0).strip() if m and "&" not in m.group(0) else None)
    block = []
    for name in ("nu_rhop_am_min", "nu_rhop_am_max"):
        block.append(re.search(rf"real\(dl\), parameter :: {name} = [^\n]*", sa).group(0))
    block.append("real(dl), parameter :: nu_rhop_fit_c = sqrt(nu_rhop_am_min*nu_rhop_am_max)")
    block.append("real(dl), parameter :: nu_rhop_fit_inv_zmax = (nu_rhop_am_max + nu_rhop_fit_c)/(nu_rhop_am_max - nu_rhop_fit_c)")
    block.append("real(dl), parameter :: nu_fit_p_denom_scale = fermi_dirac_const/((900._dl/120._dl)*zeta5)")
    m = re.search(r"real\(dl\), parameter :: nu_fit_p_c\(0:6\) = \[.*?\]", sa, re.S)
    block.append(m.group(0))
    m = re.search(r"real\(dl\), parameter :: nu_low_p_c\(0:1\) = \[.*?\]", sa, re.S)
    block.append(m.group(0))
    decl = "    ! t10a nufix: A0 pressure-fit constants\n" + "".join(f"    {x}\n" for x in block)
    anchor = "    real(dl), parameter :: am_maxp=am_max*0.9_dl\n"
    assert anchor in sb
    sb = sb.replace(anchor, anchor + decl)
    # replace the three B1 derivative functions with the A0 ones (A0 routine + wrappers)
    i0 = sb.index("    ! EFTCAMB MOD START: compatibility with massive neutrinos\n    function ThermalNuBackground_pidot(")
    i1 = sb.index("    end function ThermalNuBackground_pidotdotdot") + len("    end function ThermalNuBackground_pidotdotdot")
    j0 = sa.index("    ! EFTCAMB MOD START: compatibility with massive neutrinos\n    pure subroutine nu_pressure_lnam_derivs(")
    j1 = sa.index("    end function ThermalNuBackground_pidotdotdot") + len("    end function ThermalNuBackground_pidotdotdot")
    new = sa[j0:j1].replace("pure function", "function").replace("class(TThermalNuBackground), intent(in) :: this", "class(TThermalNuBackground) :: this")
    sb = sb[:i0] + new + sb[i1:]
    fb.write_text(sb)
    print(f"patched {fb}")


if __name__ == "__main__":
    main()
