"""Independent check of eft-par-cache-realloc-abort: different API paths, ini-based model, and pars reuse."""
import os
import subprocess
import sys

REPO = "/Users/marco/IW/eftcamb_update/EFTCAMB"
TREES = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees"


def child(tree, mode):
    os.environ["OMP_NUM_THREADS"] = "1"
    sys.path.insert(0, REPO + "/fortran/eftcamb_test/audit/checks")
    import harness
    sys.path.insert(0, tree)
    import camb
    with harness.slot(1):
        ini = REPO + "/fortran/eftcamb_test/parameters/" + mode.split(":")[1]
        pars = camb.read_ini(ini, no_validate=True)
        if "B1" in tree:
            print("    B1 tree", flush=True)
        print("    reading done", flush=True)
        r = camb.get_background(pars)
        print("    bg done", flush=True)
        h1 = r.hubble_parameter(1.0)
        print(f"    first: EFTflag={pars.EFTCAMB.EFTflag if hasattr(pars, 'EFTCAMB') else '?'} H(z=1)={h1:.12f}", flush=True)
        kind = mode.split(":")[0]
        if kind == "same_pars":
            r2 = camb.get_background(pars)
        elif kind == "calc_bg_copy":
            r2 = camb.CAMBdata()
            r2.calc_background(r.Params.copy())
        elif kind == "get_results_copy":
            p2 = r.Params.copy()
            p2.WantTransfer = False
            r2 = camb.get_background(p2)
        print(f"    second ({kind}): H(z=1)={r2.hubble_parameter(1.0):.12f}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "_child":
        child(sys.argv[2], sys.argv[3])
        sys.exit(0)
    for tree in sys.argv[1:]:
        for ini in ("2_PEFT_Omega_lin_1.ini", "4_Des_F_R_LCDM_1.ini"):
            for kind in ("same_pars", "calc_bg_copy", "get_results_copy"):
                p = subprocess.run([sys.executable, __file__, "_child", TREES + "/" + tree, f"{kind}:{ini}"],
                                   capture_output=True, text=True, timeout=600)
                print(f"{tree} {ini} {kind}: exit {p.returncode}")
                print(p.stdout.rstrip())
                for ln in p.stderr.splitlines():
                    if "rror" in ln or ln.startswith("At line"):
                        print("    " + ln)
