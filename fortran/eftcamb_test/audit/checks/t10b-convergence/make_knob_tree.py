"""Build the t10b 'knob' tree: A0 + environment-gated experiments (inert when the variables are unset).

AUDIT_T10B_PISTEP=1   : limit the scalar RK45 step to 1.5/omega_pi/IntTolBoost while EFT is active, with
                        omega_pi = sqrt(|(C + k^2 D1 + k^4 D2)/(A1 + k^2 A2)|) + |(B1 + k^2 B2)/(A1 + k^2 A2)|
                        from the last derivs call (equations.f90 GaugeInterface_ScalEv).
AUDIT_T10B_QMIN=f     : multiply the minimum source wavenumber qmin (cmbmain InitVars) by f.
AUDIT_T10B_DTAU0=f    : multiply the late-time source time step dtau0 (results SetTimeSteps) by f.
AUDIT_T10B_TON=R      : add CMB source time steps that bracket the EFT perturbation turn-on time tightly
                        (tau_on*(1 -+ 1e-7)) and a dense grid (dtau0/10) over the following R Mpc (R<=1.5: 300 Mpc).
AUDIT_T10B_PRE=R      : with AUDIT_T10B_TON, also a dense grid (dtau0/10) over the R Mpc before the turn-on.
AUDIT_T10B_EXACTON=1  : always integrate up to the EFT turn-on time before switching (no smallTime shortcut),
                        so the EFT initial conditions are never applied before EFTCAMB_turn_on_time.
AUDIT_T10B_ICREG=eps  : EFT_IC_type=1: use pi = pi' = 0 when |C + k^2 D1 + k^4 D2| < eps (|C| + k^2|D1| + k^4|D2|)
                        (regularises the attractor IC near the zero of the pi mass term, 09_EFTCAMB_IC.f90).
AUDIT_T10B_KEEPMAX=1  : RungeKuttaDP45 ind=1 keeps the caller's max_step_size instead of resetting it to 20 Mpc.
AUDIT_T10B_ICTIME=1   : start every scalar mode at min(GetTauStart(q), 0.9 tau(EFT turn-on)) so that the EFT
                        initial conditions are always applied at the same (k-independent) time, the turn-on.
python make_knob_tree.py DEST
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ENVREAD = '''
    function t10b_env(name, default) result(val)
    character(len=*), intent(in) :: name
    real(dl), intent(in) :: default
    real(dl) :: val
    character(len=64) :: s
    integer :: st, ios
    val = default
    call get_environment_variable(name, s, status=st)
    if (st == 0) then
        read(s, *, iostat=ios) val
        if (ios /= 0) val = default
    end if
    end function t10b_env
'''


def sub(text: str, old: str, new: str, path: str) -> str:
    if text.count(old) != 1:
        raise SystemExit(f"patch anchor not unique/found in {path}: {old[:60]!r}")
    return text.replace(old, new)


def main() -> None:
    dest = Path(sys.argv[1])
    if not dest.exists():
        import shutil
        shutil.copytree(harness.tree_path("A0"), dest, symlinks=True)
        (dest / "camb" / "camblib.so").unlink(missing_ok=True)
    f = dest / "fortran" / "equations.f90"
    s = f.read_text()
    s = sub(s, '''    subroutine GaugeInterface_ScalEv(EV, y, tau, tauend, tol1, ind, rk_settings, w)
    type(EvolutionVars) EV
    type(RungeKuttaDP45Settings), intent(inout) :: rk_settings
    real(dl) w(EV%nvar, 9), y(EV%nvar), tol1, tau, tauend
    integer ind
''', '''    subroutine GaugeInterface_ScalEv(EV, y, tau, tauend, tol1, ind, rk_settings, w)
    type(EvolutionVars) EV
    type(RungeKuttaDP45Settings), intent(inout) :: rk_settings
    real(dl) w(EV%nvar, 9), y(EV%nvar), tol1, tau, tauend
    integer ind
    real(dl) :: pi_A, pi_omega
    real(dl), save :: t10b_pistep = -1._dl
''', str(f))
    s = sub(s, '''        rk_settings%max_step_size = 20._dl/CP%Accuracy%IntTolBoost
    end if
    rk_settings%reuse_first_derivative = .true.

    call RungeKuttaDP45(EV, EV%ScalEqsToPropagate, derivs''', '''        rk_settings%max_step_size = 20._dl/CP%Accuracy%IntTolBoost
    end if
    if (t10b_pistep < 0._dl) t10b_pistep = t10b_env('AUDIT_T10B_PISTEP', 0._dl)
    if (t10b_pistep > 0._dl .and. CP%EFTCAMB%EFTflag /= 0 .and. EV%EFTCAMBactive) then
        pi_A = EV%eft_cache%EFTpiA1 + EV%k2_buf*EV%eft_cache%EFTpiA2
        if (pi_A /= 0._dl) then
            pi_omega = sqrt(abs((EV%eft_cache%EFTpiC + EV%k2_buf*EV%eft_cache%EFTpiD1 + &
                EV%k2_buf**2*EV%eft_cache%EFTpiD2)/pi_A)) + &
                abs((EV%eft_cache%EFTpiB1 + EV%k2_buf*EV%eft_cache%EFTpiB2)/pi_A)
            if (pi_omega > 0._dl) rk_settings%max_step_size = min(rk_settings%max_step_size, &
                1.5_dl/pi_omega/CP%Accuracy%IntTolBoost)
        end if
    end if
    rk_settings%reuse_first_derivative = .true.

    call RungeKuttaDP45(EV, EV%ScalEqsToPropagate, derivs''', str(f))
    s = sub(s, '''    end subroutine GaugeInterface_ScalEv
''', '''    end subroutine GaugeInterface_ScalEv
''' + ENVREAD, str(f))
    s = sub(s, '''    real(dl) noSwitch, smallTime, switchBoost, latePhotSwitchBoost
''', '''    real(dl) noSwitch, smallTime, switchBoost, latePhotSwitchBoost
    real(dl), save :: t10b_exacton = -1._dl
''', str(f))
    s = sub(s, '''        if (next_switch > tau + smallTime) then
''', '''        if (t10b_exacton < 0._dl) t10b_exacton = t10b_env('AUDIT_T10B_EXACTON', 0._dl)
        if (next_switch > tau + smallTime .or. (t10b_exacton > 0._dl .and. &
            next_switch == EV%EFTturnOnTime .and. next_switch > tau)) then
''', str(f))
    f.write_text(s)

    f = dest / "fortran" / "eftcamb" / "09_EFTCAMB_IC.f90"
    s = f.read_text()
    s = sub(s, """    real(dl) :: Omega0, Gamma10, Gamma20, Gamma30
""", """    real(dl) :: Omega0, Gamma10, Gamma20, Gamma30
    real(dl) :: t10b_den, t10b_scale, t10b_eps
    character(len=64) :: t10b_s
    integer :: t10b_st
""", str(f))
    s = sub(s, """        case (1)
            if (EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2 /= 0._dl) then
""", """        case (1)
            t10b_den = EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2
            t10b_scale = abs(EV%eft_cache%EFTpiC) + k2*abs(EV%eft_cache%EFTpiD1) + k2**2*abs(EV%eft_cache%EFTpiD2)
            t10b_eps = 0._dl
            call get_environment_variable('AUDIT_T10B_ICREG', t10b_s, status=t10b_st)
            if (t10b_st == 0) read(t10b_s, *) t10b_eps
            if (t10b_eps > 0._dl .and. abs(t10b_den) < t10b_eps*t10b_scale) then
                y(EV%w_ix) = 0._dl
                y(EV%w_ix + 1) = 0._dl
            else if (EV%eft_cache%EFTpiC + k2*EV%eft_cache%EFTpiD1 + k2**2*EV%eft_cache%EFTpiD2 /= 0._dl) then
""", str(f))
    f.write_text(s)

    f = dest / "fortran" / "RungeKuttaDP45.f90"
    s = f.read_text()
    s = sub(s, """    logical :: resume_after_interrupt1, resume_after_interrupt2
    real :: EV
""", """    logical :: resume_after_interrupt1, resume_after_interrupt2
    real :: EV
    real(dl) :: t10b_keep
    real(dl), save :: t10b_keepmax = -1._dl
    character(len=16) :: t10b_s
    integer :: t10b_st
""", str(f))
    s = sub(s, """            settings = RungeKuttaDP45Settings(reuse_first_derivative=settings%reuse_first_derivative)
""", """            t10b_keep = settings%max_step_size
            settings = RungeKuttaDP45Settings(reuse_first_derivative=settings%reuse_first_derivative)
            if (t10b_keepmax < 0._dl) then
                t10b_keepmax = 0._dl
                call get_environment_variable('AUDIT_T10B_KEEPMAX', t10b_s, status=t10b_st)
                if (t10b_st == 0) t10b_keepmax = 1._dl
            end if
            if (t10b_keepmax > 0._dl) settings%max_step_size = t10b_keep
""", str(f))
    f.write_text(s)

    f = dest / "fortran" / "cmbmain.f90"
    s = f.read_text()
    s = sub(s, '''    real(dl) taumin, maxq, initAccuracyBoost
    integer itf
''', '''    real(dl) taumin, maxq, initAccuracyBoost
    integer itf
    real(dl) :: t10b_f
    character(len=64) :: t10b_s
    integer :: t10b_st
''', str(f))
    s = sub(s, '''        qmin = qmin0/State%tau0/initAccuracyBoost
    else''', '''        qmin = qmin0/State%tau0/initAccuracyBoost
        call get_environment_variable('AUDIT_T10B_QMIN', t10b_s, status=t10b_st)
        if (t10b_st == 0) then
            read(t10b_s, *) t10b_f
            qmin = qmin*t10b_f
        end if
    else''', str(f))
    s = sub(s, '''    w = 0
    y = 0
    call initial(EV, y, taustart)
    if (global_error_flag /= 0) return

    tau = taustart
    ind = 1
''', '''    w = 0
    y = 0
    t10b_tau0 = taustart
    call get_environment_variable('AUDIT_T10B_ICTIME', t10b_s, status=t10b_st)
    if (t10b_st == 0 .and. CP%EFTCAMB%EFTflag /= 0) then
        t10b_ton = State%DeltaTime(0._dl, CP%EFTCAMB%EFTCAMB_pert_turn_on)
        if (0.9_dl*t10b_ton >= State%ThermoData%tauminn) t10b_tau0 = min(taustart, 0.9_dl*t10b_ton)
    end if
    call initial(EV, y, t10b_tau0)
    if (global_error_flag /= 0) return

    tau = t10b_tau0
    ind = 1
''', str(f))
    s = sub(s, '''    real(dl) w(EV%nvar, 9), y(EV%nvar), sources(ThisSources%SourceNum)
''', '''    real(dl) w(EV%nvar, 9), y(EV%nvar), sources(ThisSources%SourceNum)
    character(len=64) :: t10b_s
    integer :: t10b_st
    real(dl) :: t10b_ton, t10b_tau0
''', str(f))
    f.write_text(s)

    f = dest / "fortran" / "results.f90"
    s = f.read_text()
    s = sub(s, '''    integer nri0, nstep
''', '''    integer nri0, nstep
    real(dl) :: t10b_f, t10b_r
    character(len=64) :: t10b_s
    integer :: t10b_st
''', str(f))
    s = sub(s, '''    call TimeSteps%Add_delta(State%taurend, State%tau0, dtau0)
''', '''    call get_environment_variable('AUDIT_T10B_DTAU0', t10b_s, status=t10b_st)
    if (t10b_st == 0) then
        read(t10b_s, *) t10b_f
        dtau0 = dtau0*t10b_f
    end if
    call TimeSteps%Add_delta(State%taurend, State%tau0, dtau0)
    call get_environment_variable('AUDIT_T10B_TON', t10b_s, status=t10b_st)
    if (t10b_st == 0 .and. State%CP%EFTCAMB%EFTflag /= 0) then
        read(t10b_s, *) t10b_r
        if (t10b_r <= 1.5_dl) t10b_r = 300._dl
        t10b_f = State%DeltaTime(0._dl, State%CP%EFTCAMB%EFTCAMB_pert_turn_on)
        if (t10b_f > State%taurst .and. t10b_f < State%tau0) then
            call TimeSteps%Add(t10b_f*(1 - 1.e-7_dl), t10b_f*(1 + 1.e-7_dl), 1)
            call TimeSteps%Add_delta(t10b_f*(1 + 1.e-7_dl), min(State%tau0, t10b_f + t10b_r), dtau0/10)
            call get_environment_variable('AUDIT_T10B_PRE', t10b_s, status=t10b_st)
            if (t10b_st == 0) then
                read(t10b_s, *) t10b_r
                call TimeSteps%Add_delta(max(State%taurend, t10b_f - t10b_r), t10b_f*(1 - 1.e-7_dl), dtau0/10)
            end if
        end if
    end if
''', str(f))
    f.write_text(s)
    print("patched", dest)
    if "--build" in sys.argv:
        print(harness.build_tree(dest, dest))


if __name__ == "__main__":
    main()
