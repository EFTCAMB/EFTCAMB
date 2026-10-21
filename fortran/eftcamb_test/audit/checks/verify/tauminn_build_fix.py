"""Verifier (t3-eft-turnon-before-tauminn-lmax): build a copy of A0 with two env-var-selected candidate behaviours.

AUDIT_TAUMIN_MODE unset/0: original code.
AUDIT_TAUMIN_MODE=1: proposed fix - InitVars starts the thermal history at min(GetTauStart(maxq), tau(EFT turn-on)).
AUDIT_TAUMIN_MODE=2: guard removed - equations.f90 initial() does not raise the 'starts before thermo tauminn' error.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "tauminn_fix"


def patch(path: Path, old: str, new: str) -> None:
    s = path.read_text()
    if new in s:
        return
    assert s.count(old) == 1, (path, old)
    path.write_text(s.replace(old, new))


def main() -> None:
    if not DEST.exists():
        import shutil

        shutil.copytree(harness.TREES["A0"], DEST, symlinks=True)
    patch(DEST / "fortran" / "cmbmain.f90", """    taumin = GetTauStart(maxq)
""", """    taumin = GetTauStart(maxq)
    ! AUDIT candidate fix (AUDIT_TAUMIN_MODE=1)
    block
        character(len=16) :: audit_env
        integer :: audit_st
        real(dl) :: audit_tau_on
        call get_environment_variable('AUDIT_TAUMIN_MODE', audit_env, status=audit_st)
        if (audit_st == 0 .and. trim(audit_env) == '1' .and. CP%EFTCAMB%EFTflag /= 0) then
            audit_tau_on = State%DeltaTime(0._dl, CP%EFTCAMB%EFTCAMB_pert_turn_on)
            write(*, *) 'AUDIT InitVars taumin, tau_on = ', taumin, audit_tau_on
            if (audit_tau_on > 0._dl) taumin = min(taumin, audit_tau_on)
        end if
    end block
""")
    patch(DEST / "fortran" / "equations.f90", """    if (CP%EFTCAMB%EFTflag /= 0 .and. EV%EFTturnOnTime < EV%ThermoData%tauminn) then
        call GlobalError(FormatString('EFTCAMB starts before thermo tauminn, ' // &
            'EFT_pert_turn_on, EFTturnOnTime, tauminn = %f, %f, %f.', CP%EFTCAMB%EFTCAMB_pert_turn_on, &
            EV%EFTturnOnTime, EV%ThermoData%tauminn))
        return
    end if
""", """    block
        character(len=16) :: audit_env
        integer :: audit_st
        call get_environment_variable('AUDIT_TAUMIN_MODE', audit_env, status=audit_st)
        if (.not. (audit_st == 0 .and. trim(audit_env) == '2')) then
            if (CP%EFTCAMB%EFTflag /= 0 .and. EV%EFTturnOnTime < EV%ThermoData%tauminn) then
                call GlobalError(FormatString('EFTCAMB starts before thermo tauminn, ' // &
                    'EFT_pert_turn_on, EFTturnOnTime, tauminn = %f, %f, %f.', CP%EFTCAMB%EFTCAMB_pert_turn_on, &
                    EV%EFTturnOnTime, EV%ThermoData%tauminn))
                return
            end if
        end if
    end block
""")
    print(harness.build_tree(DEST, DEST))


if __name__ == "__main__":
    main()
