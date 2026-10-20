"""Attribution test for finding evolve-metric-h-inconsistent (F-065): instrumented tree + exact massive-neutrino delta p.

Copies work/verify/mh_instr (built by metric_h_build_instr.py) to work/verify/mh_instr_dgp and replaces the EFTCAMB
adiabatic approximation in MassiveNuVars (dgp = gpnu_t*clxnu, i.e. delta p = w delta rho) by the exact pressure
perturbation: Nu_Integrate_L012 (dpnu, already divided by 3 there -> delta p in massless-nu units) for the full
hierarchy, y(nu_ix+1) = delta p/rho for the approximate hierarchy (see SwitchToMassiveNuApprox and the approx density
equation, equations.f90:3576). Builds it (retries for the forutils parallel-make race).

Usage: python metric_h_build_dgp.py
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

SRC = harness.AUDIT_ROOT / "work" / "verify" / "mh_instr"
DEST = harness.AUDIT_ROOT / "work" / "verify" / "mh_instr_dgp"


def patch(path: Path, old: str, new: str) -> None:
    s = path.read_text()
    assert s.count(old) == 1, (path, old[:60], s.count(old))
    path.write_text(s.replace(old, new))


if not DEST.exists():
    shutil.copytree(SRC, DEST, symlinks=True)
    eq = DEST / "fortran" / "equations.f90"
    patch(eq, "    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu\n\n    ! EFTCAMB MOD START: compatibility",
          "    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu\n    real(dl) dpnu_x, pinu_x\n\n"
          "    ! EFTCAMB MOD START: compatibility")
    patch(eq,
          "        if (EV%MassiveNuApprox(nu_i)) then\n"
          "            clxnu = y(EV%nu_ix(nu_i))\n"
          "            qnu = y(EV%nu_ix(nu_i) + 2)\n"
          "        else\n"
          "            ! Integrate over q\n"
          "            call Nu_Integrate_L01(EV, y, a, nu_i, clxnu, qnu)\n",
          "        dpnu_x = 0._dl\n"
          "        if (EV%MassiveNuApprox(nu_i)) then\n"
          "            clxnu = y(EV%nu_ix(nu_i))\n"
          "            qnu = y(EV%nu_ix(nu_i) + 2)\n"
          "            dpnu_x = y(EV%nu_ix(nu_i) + 1)\n"
          "        else\n"
          "            ! Integrate over q\n"
          "            if (present(dgp)) then\n"
          "                call Nu_Integrate_L012(EV, y, a, nu_i, clxnu, qnu, dpnu_x, pinu_x)\n"
          "                dpnu_x = dpnu_x/rhonu\n"
          "            else\n"
          "                call Nu_Integrate_L01(EV, y, a, nu_i, clxnu, qnu)\n"
          "            end if\n")
    patch(eq, "            dgp = dgp + gpnu_t*clxnu\n", "            dgp = dgp + grhonu_t*dpnu_x\n")
    print("patched", DEST)
for attempt in range(4):
    r = harness.build_tree(DEST, DEST)
    log = Path(r["log"]).read_text(errors="replace")
    ok = "Error" not in log and "camblib.so ../camb/" in log
    print(attempt, r["elapsed"], "clean" if ok else "failed")
    if ok:
        break
