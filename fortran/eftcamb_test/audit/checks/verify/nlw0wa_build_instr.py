"""Verify nonlinear-halofit-w0wa-lcdm-for-eft: build an instrumented copy of A0.

halofit.f90 in the copy prints the (w, wa) passed to Halofit and HMcode (to stderr, prefixed "AUDIT") and lets the
environment variables EFT_AUDIT_HF_W / EFT_AUDIT_HF_WA override them. Nothing else changes.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
DEST = SCR / "work/verify/nlw0wa/instr"
SRC = SCR / "trees/A0"


def block(wvar, wavar, tag):
    return (
        "                block\n"
        "                    character(len=64) :: audit_buf\n"
        "                    integer :: audit_st\n"
        "                    call get_environment_variable('EFT_AUDIT_HF_W', audit_buf, status=audit_st)\n"
        f"                    if (audit_st == 0) read(audit_buf, *) {wvar}\n"
        "                    call get_environment_variable('EFT_AUDIT_HF_WA', audit_buf, status=audit_st)\n"
        f"                    if (audit_st == 0) read(audit_buf, *) {wavar}\n"
        f"                    write(0, '(a,2es16.8)') 'AUDIT {tag} w wa =', {wvar}, {wavar}\n"
        "                end block\n"
    )


def patch(path: Path):
    s = path.read_text()
    a = ("                ! Original code:\n                ! call Params%DarkEnergy%Effective_w_wa(w_eff, wa_eff)\n"
         "                ! EFTCAMB MOD END\n")
    assert s.count(a) == 1, "halofit anchor"
    s = s.replace(a, a + block("w_eff", "wa_eff", "halofit"))
    b = ("        ! Original code:\n        ! call CP%DarkEnergy%Effective_w_wa(cosm%w, cosm%wa)\n"
         "        ! EFTCAMB MOD END\n")
    assert s.count(b) == 1, "hmcode anchor"
    s = s.replace(b, b + block("cosm%w", "cosm%wa", "hmcode"))
    path.write_text(s)


if __name__ == "__main__":
    if not DEST.exists():
        import shutil

        shutil.copytree(SRC, DEST, symlinks=True)
        patch(DEST / "fortran/halofit.f90")
    print(harness.build_tree(SRC, DEST))
