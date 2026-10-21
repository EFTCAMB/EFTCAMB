"""Attribution of curvature-silently-wrong-with-eft with a DIAGNOSTIC patch of a copy of A0 (not a proposed fix).

Different code path from the earlier t3 verification (which patched adotoa inside derivs): here the curvature is put
back at MODEL level, i.e. what the finding's alternative fix proposes, and the sources are switched separately.
Runtime switches, env AUDIT_CURVSIL (bitmask, read once per routine):
  bit 1: Pure-EFT designer (007p1) compute_adotoa adds eft_par_cache%grhok (H^2 = (grho + grhok)/3) and
         compute_H_derivs adds grhok/6 to Hdot (true dH/dtau in a curved FRW)
  bit 2: curved cothxor in derivs' EFT branch (equations.f90 ~2907) instead of 1/tau
  bit 4: in the EFTflag/=0, EFT-inactive temperature source use CAMB's ISW = 2*phidot*exptau instead of EFTCAMB's
         flat-space rewrite -2 sigma (Hdot - 2 H^2) - 2 H etak + dgq - diff_rhopi/k (equations.f90 ~3876)
AUDIT_CURVSIL=0 must reproduce A0.
Usage: python curvsil_attrib.py build | run | analyze
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify"
DEST = WORK / "curvsil_tree"
OUT = WORK / "curvsil" / "attrib"
REF = WORK / "curvsil"
P = harness.PARAMS_DIR
COMMON = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "2200"}

READ_ENV = """    if (audit_flag < 0) then
        call get_environment_variable('AUDIT_CURVSIL', audit_env)
        audit_flag = 0
        if (len_trim(audit_env) > 0) read(audit_env, *) audit_flag
    end if
"""
DECL = "    integer, save :: audit_flag = -1\n    character(len=32) :: audit_env\n"

EDITS_007 = [
    (
        "    subroutine EFTCAMBPureEFTstdComputeAdotoa(self, a, eft_par_cache, eft_cache)\n\n    implicit none\n",
        "    subroutine EFTCAMBPureEFTstdComputeAdotoa(self, a, eft_par_cache, eft_cache)\n\n    implicit none\n" + DECL,
    ),
    (
        "    eft_cache%adotoa = sqrt((eft_cache%grhom_t + eft_cache%grhov_t)/3._dl)\n\n"
        "    end subroutine EFTCAMBPureEFTstdComputeAdotoa\n",
        "    eft_cache%adotoa = sqrt((eft_cache%grhom_t + eft_cache%grhov_t)/3._dl)\n" + READ_ENV +
        "    if (iand(audit_flag, 1) /= 0) eft_cache%adotoa = sqrt((eft_cache%grhom_t + eft_cache%grhov_t + "
        "eft_par_cache%grhok)/3._dl) ! AUDIT CURVSIL\n\n"
        "    end subroutine EFTCAMBPureEFTstdComputeAdotoa\n",
    ),
    (
        "    subroutine EFTCAMBPureEFTstdComputeHubbleDer(self, a, eft_par_cache, eft_cache)\n\n    implicit none\n",
        "    subroutine EFTCAMBPureEFTstdComputeHubbleDer(self, a, eft_par_cache, eft_cache)\n\n    implicit none\n" + DECL,
    ),
    (
        "    eft_cache%Hdot = -0.5_dl*(eft_cache%adotoa**2 + eft_cache%gpresm_t + eft_cache%gpiv_t)\n",
        "    eft_cache%Hdot = -0.5_dl*(eft_cache%adotoa**2 + eft_cache%gpresm_t + eft_cache%gpiv_t)\n" + READ_ENV +
        "    if (iand(audit_flag, 1) /= 0) eft_cache%Hdot = eft_cache%Hdot + eft_par_cache%grhok/6._dl ! AUDIT CURVSIL\n",
    ),
]

EDITS_EQ = [
    ("    real(dl) cothxor ! 1/tau in flat case\n", "    real(dl) cothxor ! 1/tau in flat case\n" + DECL),
    ("    ! EFTCAMB MOD START: compute adotoa\n", READ_ENV + "    ! EFTCAMB MOD START: compute adotoa\n"),
    (
        "        ! call CP%EFTCAMB%model%EFTCAMBModelComputeQSMuSigma( a, k, CP%eft_par_cache , EV%eft_cache ) minsup\n"
        "        !\n"
        "        cothxor = 1._dl/tau\n",
        "        ! call CP%EFTCAMB%model%EFTCAMBModelComputeQSMuSigma( a, k, CP%eft_par_cache , EV%eft_cache ) minsup\n"
        "        !\n"
        "        cothxor = 1._dl/tau\n"
        "        if (.not. State%flat .and. iand(audit_flag, 2) /= 0) & ! AUDIT CURVSIL\n"
        "            cothxor = 1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius\n",
    ),
    (
        "                if (CP%EFTCAMB%EFTflag /= 0) then\n"
        "                    ISW = -2._dl*sigma*(EV%eft_cache%Hdot - 2._dl*adotoa**2) - 2._dl*adotoa*etak + dgq - "
        "diff_rhopi/k\n",
        "                if (CP%EFTCAMB%EFTflag /= 0 .and. iand(audit_flag, 4) == 0) then ! AUDIT CURVSIL\n"
        "                    ISW = -2._dl*sigma*(EV%eft_cache%Hdot - 2._dl*adotoa**2) - 2._dl*adotoa*etak + dgq - "
        "diff_rhopi/k\n",
    ),
]


def patch(path: Path, edits) -> None:
    s = path.read_text()
    if "AUDIT CURVSIL" in s:
        return
    for old, new in edits:
        assert s.count(old) == 1, (path.name, old[:80], s.count(old))
        s = s.replace(old, new)
    path.write_text(s)


def build() -> None:
    import shutil

    if not DEST.exists():
        shutil.copytree(harness.tree_path("A0"), DEST, symlinks=True)
    patch(DEST / "fortran" / "eftcamb" / "07f_designer_models" / "007p1_Pure_EFT_std.f90", EDITS_007)
    patch(DEST / "fortran" / "equations.f90", EDITS_EQ)
    so = DEST / "camb" / "camblib.so"
    so.unlink(missing_ok=True)  # the copied A0 library would make a failed build look ok
    for _ in range(4):  # forutils' Makefile races under make -j4 on a fresh build dir: retry
        r = harness.build_tree(DEST, DEST)
        print(r)
        if r["ok"]:
            break


def cases():
    c = []
    for omk in ("0.02", "-0.02"):
        for flag in (0, 1, 2, 3, 4, 5, 7):
            c.append((P / "2_PEFT_wDE_LCDM_1.ini", f"PEFT_omk{omk}_f{flag}", omk, flag))
    c.append((P / "2_PEFT_wDE_LCDM_1.ini", "PEFT_omk0_f7", "0", 7))  # patch inert when flat
    c.append((P / "1_EFT_GR.ini", "GR_omk0.02_f7", "0.02", 7))  # patch inert for EFTflag=0
    return c


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    def one(cs):
        ini, name, omk, flag = cs
        if (OUT / f"{name}.npz").exists():
            return name, "cached"
        r = harness.run_model(DEST, ini, OUT / f"{name}.npz", overrides=dict(COMMON, omk=omk), threads=1,
                              timeout=900, extra_env={"AUDIT_CURVSIL": str(flag)})
        return name, r["status"]

    with ThreadPoolExecutor(4) as ex:
        for name, st in ex.map(one, cases()):
            print(name, st, flush=True)


def cmp(a, b):
    ca, cb = a["cl_total"], b["cl_total"]
    n = min(len(ca), len(cb), 2201)
    r = ca[2:n] / cb[2:n] - 1
    te = (ca[2:n, 3] - cb[2:n, 3]) / np.sqrt(cb[2:n, 0] * cb[2:n, 1])
    return {
        "TT_l2": float(r[0, 0]), "TT_l10": float(r[8, 0]), "TT_max_l>=30": float(np.max(np.abs(r[28:, 0]))),
        "EE_max": float(np.max(np.abs(r[:, 1]))), "TE_max_norm": float(np.max(np.abs(te))),
        "sigma8": float(a["sigma8"][-1] / b["sigma8"][-1] - 1),
        "pk_z0_max": float(np.max(np.abs(a["pk"][-1] / b["pk"][-1] - 1))),
    }


def analyze() -> None:
    rows = {}
    for ini, name, omk, flag in cases():
        f = OUT / f"{name}.npz"
        if not f.exists():
            continue
        d = harness.load(f)
        ref = harness.load(REF / f"B0_GR_omk{omk}.npz")
        rows[f"{name} vs B0 CAMB"] = cmp(d, ref)
        if flag == 0:
            rows[f"{name} vs frozen A0 (must be 0)"] = cmp(d, harness.load(REF / f"A0_PEFT_omk{omk}.npz"))
        if name.startswith("GR_"):
            rows[f"{name} vs frozen A0 GR (must be 0)"] = cmp(d, harness.load(REF / f"A0_GR_omk{omk}.npz"))
    print(f"{'case':38s} {'TT_l2':>10s} {'TT_l10':>10s} {'TTl>=30':>9s} {'EE_max':>9s} {'TEnorm':>9s} "
          f"{'sigma8':>10s} {'pk_z0':>9s}")
    for k, v in rows.items():
        print(f"{k:38s} {v['TT_l2']:+10.3e} {v['TT_l10']:+10.3e} {v['TT_max_l>=30']:9.2e} {v['EE_max']:9.2e} "
              f"{v['TE_max_norm']:9.2e} {v['sigma8']:+10.3e} {v['pk_z0_max']:9.2e}")
    (OUT / "curvsil_attrib_summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    {"build": build, "run": run, "analyze": analyze}[sys.argv[1]]()
