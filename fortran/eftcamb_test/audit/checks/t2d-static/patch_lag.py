"""Make the LAG tree copy: FIX3 (A0 + Extended Galileon neutrino/bisection fix) plus instrumentation in derivs.

For full-mapping models derivs calls compute_background_EFT_functions BEFORE compute_adotoa / compute_H_derivs
(equations.f90:2813-2817) and never again, so model background EFT functions that read eft_cache%adotoa, Hdot or
Hdotdot (Extended Galileon, Beyond Horndeski) see the values left by the previous derivs call of that k mode.
The instrumentation, after the second compute_H_derivs call (full models only):
  * saves the cache, re-evaluates compute_background_EFT_functions with the now-current adotoa/Hdot/Hdotdot,
  * writes to stderr a line 'LAGMAX field a k tau used fresh rel' whenever the relative difference between the value
    used by derivs (lagged) and the fresh value reaches a new maximum (a >= 1e-3), and a 'LAGSUM' line every 20000 calls,
  * keeps the fresh values when the environment variable EFT_LAG_FIX=1 (candidate fix), otherwise restores the saved
    cache (bit-identical to FIX3 behaviour).
Run with OMP_NUM_THREADS=1 (module-level counters, not thread safe by design).

Usage: python patch_lag.py [--build]
"""

from __future__ import annotations

import difflib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

SRC = harness.AUDIT_ROOT / "work/t2d-static/FIX3"
DEST = harness.AUDIT_ROOT / "work/t2d-static/LAG"
REL = "fortran/equations.f90"

EDITS = [
    ("    integer, parameter :: basic_num_eqns = 4\n",
     "    integer, parameter :: basic_num_eqns = 4\n"
     "    ! AUDIT LAG instrumentation (module state, single thread only)\n"
     "    real(dl), save :: lag_max(4) = 0._dl\n"
     "    integer, save :: lag_calls = 0\n"
     "    logical, save :: lag_init = .false., lag_fix = .false.\n"),
    ("    real(dl) cothxor ! 1/tau in flat case\n    real(dl) xe, Trad, Delta_TM, Tmat, Delta_TCMB\n",
     "    real(dl) cothxor ! 1/tau in flat case\n    real(dl) xe, Trad, Delta_TM, Tmat, Delta_TCMB\n"
     "    real(dl) lag_used(4), lag_fresh(4), lag_rel(4)\n    integer lag_i\n    character(len=8) lag_env\n"
     "    type(TEFTCAMB_timestep_cache) :: lag_cache\n"),
    ("        ! compute third H derivate:\n        call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)\n"
     "        ! store:\n        adotdota = EV%eft_cache%Hdot + EV%eft_cache%adotoa**2\n",
     "        ! compute third H derivate:\n        call CP%EFTCAMB%model%compute_H_derivs(a, CP%eft_par_cache, EV%eft_cache)\n"
     "        ! AUDIT LAG instrumentation:\n"
     "        if (.not. CP%EFTCAMB%EFTCAMB_model_is_designer) then\n"
     "            if (.not. lag_init) then\n"
     "                lag_env = ''\n"
     "                call get_environment_variable('EFT_LAG_FIX', lag_env)\n"
     "                lag_fix = trim(lag_env) == '1'\n"
     "                lag_init = .true.\n"
     "                write(0, '(a,l2)') 'LAGINIT fix =', lag_fix\n"
     "            end if\n"
     "            lag_cache = EV%eft_cache\n"
     "            lag_used = [EV%eft_cache%EFTc, EV%eft_cache%EFTLambda, EV%eft_cache%EFTcdot, EV%eft_cache%EFTLambdadot]\n"
     "            call CP%EFTCAMB%model%compute_background_EFT_functions(a, CP%eft_par_cache, EV%eft_cache)\n"
     "            lag_fresh = [EV%eft_cache%EFTc, EV%eft_cache%EFTLambda, EV%eft_cache%EFTcdot, &\n"
     "                EV%eft_cache%EFTLambdadot]\n"
     "            lag_calls = lag_calls + 1\n"
     "            do lag_i = 1, 4\n"
     "                lag_rel(lag_i) = abs(lag_used(lag_i) - lag_fresh(lag_i))/max(abs(lag_fresh(lag_i)), 1.d-300)\n"
     "                if (a >= 1.d-3 .and. lag_rel(lag_i) > lag_max(lag_i)) then\n"
     "                    lag_max(lag_i) = lag_rel(lag_i)\n"
     "                    write(0, '(a,i2,6es14.5)') 'LAGMAX ', lag_i, a, k, tau, lag_used(lag_i), lag_fresh(lag_i), &\n"
     "                        lag_rel(lag_i)\n"
     "                end if\n"
     "            end do\n"
     "            if (mod(lag_calls, 20000) == 0) write(0, '(a,i9,4es12.4)') 'LAGSUM ', lag_calls, lag_max\n"
     "            if (.not. lag_fix) EV%eft_cache = lag_cache\n"
     "        end if\n"
     "        ! store:\n        adotdota = EV%eft_cache%Hdot + EV%eft_cache%adotoa**2\n"),
]


def main() -> int:
    if not DEST.exists():
        shutil.copytree(SRC, DEST, symlinks=True)
    path = DEST / REL
    old = (SRC / REL).read_text()
    new = old
    for a, b in EDITS:
        assert new.count(a) == 1, a[:80]
        new = new.replace(a, b)
    path.write_text(new)
    sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                               f"FIX3/{REL}", f"LAG/{REL}"))
    if "--build" in sys.argv:
        print(harness.build_tree(SRC, DEST, make_args=["-j1", "CLUSTER_SAFE=0"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
