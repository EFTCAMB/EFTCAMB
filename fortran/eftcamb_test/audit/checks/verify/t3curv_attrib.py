"""Attribution of t3-eft-curvature-ignored with runtime switches (one build, several variants).

DIAGNOSTIC patch of a copy of A0 (not a proposed fix), selected by env AUDIT_CURVDIAG (bitmask, read once):
  bit 1: in derivs' EFT background branch restore curvature in adotoa (+grhok/3), use CAMB's GR adotdota and the curved
         cothxor  (the stated root cause)
  bit 2: in the EFTflag/=0-but-EFT-inactive ISW source (EFTCAMB's replacement of 2*phidot) fall back to CAMB's
         2*phidot when not flat
AUDIT_CURVDIAG=0 must reproduce A0.
Usage: python t3curv_attrib.py build | run
"""

import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import t3curv_lens_physics as L  # noqa: E402

DEST = harness.AUDIT_ROOT / "work" / "verify" / "t3curv_attrib_tree"

EDITS = [
    (
        "    real(dl) cothxor ! 1/tau in flat case\n",
        "    real(dl) cothxor ! 1/tau in flat case\n"
        "    integer, save :: audit_curv = -1\n"
        "    character(len=32) :: audit_env\n",
    ),
    (
        "    ! EFTCAMB MOD START: compute adotoa\n",
        "    if (audit_curv < 0) then\n"
        "        call get_environment_variable('AUDIT_CURVDIAG', audit_env)\n"
        "        audit_curv = 0\n"
        "        if (len_trim(audit_env) > 0) read(audit_env, *) audit_curv\n"
        "    end if\n"
        "    ! EFTCAMB MOD START: compute adotoa\n",
    ),
    (
        "        ! call CP%EFTCAMB%model%EFTCAMBModelComputeQSMuSigma( a, k, CP%eft_par_cache , EV%eft_cache ) minsup\n"
        "        !\n"
        "        cothxor = 1._dl/tau\n",
        "        ! call CP%EFTCAMB%model%EFTCAMBModelComputeQSMuSigma( a, k, CP%eft_par_cache , EV%eft_cache ) minsup\n"
        "        !\n"
        "        cothxor = 1._dl/tau\n"
        "        if (.not. State%flat .and. iand(audit_curv, 1) /= 0) then\n"
        "            adotoa = sqrt(adotoa**2 + State%grhok/3._dl)\n"
        "            adotdota = (adotoa*adotoa - gpres_noDE - EV%eft_cache%gpiv_t)/2._dl\n"
        "            cothxor = 1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius\n"
        "        end if\n",
    ),
    (
        "                if (CP%EFTCAMB%EFTflag /= 0) then\n"
        "                    ISW = -2._dl*sigma*(EV%eft_cache%Hdot - 2._dl*adotoa**2) - 2._dl*adotoa*etak + dgq - diff_rhopi/k\n",
        "                if (CP%EFTCAMB%EFTflag /= 0 .and. (State%flat .or. iand(audit_curv, 2) == 0)) then\n"
        "                    ISW = -2._dl*sigma*(EV%eft_cache%Hdot - 2._dl*adotoa**2) - 2._dl*adotoa*etak + dgq - diff_rhopi/k\n",
    ),
]


def build() -> None:
    if not DEST.exists():
        shutil.copytree(harness.tree_path("A0"), DEST, symlinks=True)
    f = DEST / "fortran" / "equations.f90"
    s = f.read_text()
    if "AUDIT_CURVDIAG" not in s:
        for old, new in EDITS:
            assert s.count(old) == 1, old
            s = s.replace(old, new)
        f.write_text(s)
    for _ in range(4):  # forutils make race on a fresh build dir under -j4: retry
        r = harness.build_tree("A0", DEST)
        if r["ok"] and "camblib.so ../camb/" in Path(r["log"]).read_text():
            break
    print(r)


def jobs():
    j = []
    for mode in ("0", "1", "2", "3"):
        for omk in ("0.01", "0.05", "-0.05"):
            j.append((mode, {**L.COMMON, "omk": omk, "PureEFTmodelOmega": "0"}, f"X{mode}_PEFT0_omk{omk}"))
    j.append(("3", {**L.COMMON, "omk": "0", "PureEFTmodelOmega": "0"}, "X3_PEFT0_omk0"))
    j.append(("3", {**L.COMMON, "omk": "0.05"}, "X3_GR_omk0.05"))
    for omk in ("0", "0.01", "0.05", "-0.05"):  # EFT equations active (Omega=1e-6 above the RGR threshold)
        j.append(("3", {**L.COMMON, "omk": omk, "PureEFTmodelOmega": "1", "EFTOmega0": "1e-6"}, f"X3_PEFTom6_omk{omk}"))
        j.append(("0", {**L.COMMON, "omk": omk, "PureEFTmodelOmega": "1", "EFTOmega0": "1e-6"}, f"X0_PEFTom6_omk{omk}"))
    return j


def run() -> None:
    def one(x):
        mode, ov, name = x
        out = L.OUT / f"{name}.npz"
        if out.exists():
            return {"status": "ok"}
        ini = L.GR if "_GR_" in name else L.PEFT
        return harness.run_model(str(DEST), ini, out, overrides=ov, threads=1, timeout=900,
                                 extra_env={"AUDIT_CURVDIAG": mode})

    with ThreadPoolExecutor(4) as ex:
        for s in ex.map(one, jobs()):
            if s["status"] != "ok":
                print("FAILED", s.get("out"), s["status"], s.get("log_tail", [])[-4:])
    keys = ("TT", "TT_l>=30", "EE", "TE", "PP", "pk", "sigma8", "TT_l2_ratio")
    print(f"{'case (vs A0 EFTflag=0, same omk)':52s} " + " ".join(f"{k[:9]:>9s}" for k in keys))
    for _, _, name in jobs():
        omk = name.split("omk")[1]
        r = L.row(name, f"A0_GR_omk{omk}")
        print(f"{name:52s} " + " ".join(f"{r[k]:9.2e}" for k in keys))
    for omk in ("0.01", "0.05", "-0.05"):
        r = harness.compare(L.OUT / f"X0_PEFT0_omk{omk}.npz", L.OUT / f"A0_PEFT0_omk{omk}.npz")
        print(f"mode0 vs unpatched A0 omk={omk}: bit_identical={r['bit_identical']} worst_cl={r['worst_cl']:.2e}")


if __name__ == "__main__":
    {"build": build, "run": run}[sys.argv[1]]()
