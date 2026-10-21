"""Make the FIX3 tree copy (candidate fix for the Extended Galileon adotoa: massive neutrinos computed locally instead
of reading the not-yet-filled eft_cache%grhonu_tot, and a real bisection tolerance instead of the integer 10**(-1)=0).

Usage: python patch_fix3.py   (copies trees/A0 to work/t2d-static/FIX3, patches it, builds it with harness.build_tree)
Prints the unified diff of the patch.
"""

from __future__ import annotations

import difflib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

SRC = harness.TREES["A0"]
DEST = harness.AUDIT_ROOT / "work/t2d-static/FIX3"
REL = "fortran/eftcamb/08f_full_models/008p7_Extended_Galileon.f90"

EDITS = [
    ("    use EFTCAMB_abstract_model_full\n",
     "    use EFTCAMB_abstract_model_full\n    use MassiveNu\n"),
    ("    real(dl) :: ATemp1, ATemp2, BTemp1, BTemp2, HorizAsyntB\n\n    Omega_phi0 = eft_par_cache%omegav\n"
     "    a2 = a*a\n\n"
     "    Omega_tot = (eft_par_cache%omegac + eft_par_cache%omegab)*a**(-3.) + &\n"
     "        (eft_par_cache%omegag + eft_par_cache%omegar)*a**(-4.) + &\n"
     "        eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2.*a2)\n",
     "    real(dl) :: ATemp1, ATemp2, BTemp1, BTemp2, HorizAsyntB\n"
     "    real(dl) :: grhonu_loc, rhonu_i, pnu_i\n\n    Omega_phi0 = eft_par_cache%omegav\n"
     "    a2 = a*a\n\n"
     "    ! massive neutrinos computed here: eft_cache%grhonu_tot is filled by the callers only after compute_adotoa\n"
     "    grhonu_loc = 0._dl\n"
     "    if (eft_par_cache%Num_Nu_Massive /= 0) then\n"
     "        do nu_i = 1, eft_par_cache%Nu_mass_eigenstates\n"
     "            rhonu_i = 0._dl\n"
     "            pnu_i = 0._dl\n"
     "            call ThermalNuBackground%rho_P(a*eft_par_cache%nu_masses(nu_i), rhonu_i, pnu_i)\n"
     "            grhonu_loc = grhonu_loc + eft_par_cache%grhormass(nu_i)/a2*rhonu_i\n"
     "        end do\n"
     "    end if\n"
     "    Omega_tot = (eft_par_cache%omegac + eft_par_cache%omegab)*a**(-3) + &\n"
     "        (eft_par_cache%omegag + eft_par_cache%omegar)*a**(-4) + &\n"
     "        grhonu_loc/(3._dl*eft_par_cache%h0_Mpc**2*a2)\n"),
    ("    do while (sqrt(fsoluction**2) > 10**(-1) .and. counter < 50**1)\n",
     "    do while (abs(fsoluction) > 1.e-13_dl .and. counter < 100)\n"),
]


def main() -> None:
    if not DEST.exists():
        shutil.copytree(SRC, DEST, symlinks=True)
    f = DEST / REL
    old = (SRC / REL).read_text()
    new = old
    for a, b in EDITS:
        assert new.count(a) == 1, a[:80]
        new = new.replace(a, b)
    f.write_text(new)
    sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), new.splitlines(True), f"A0/{REL}", f"FIX3/{REL}"))
    if "--build" in sys.argv:
        print(harness.build_tree(SRC, DEST))


if __name__ == "__main__":
    main()
