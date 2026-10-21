"""Independent check: use module-level camb.write_ini(validate=False) and inspect the raw written file
for every EFT key present in the original flattened ini; also check DarkEnergy settings written."""
import os, sys, glob
os.environ["OMP_NUM_THREADS"] = "1"
REPO = "/Users/marco/IW/eftcamb_update/EFTCAMB"
sys.path.insert(0, REPO + "/fortran/eftcamb_test/audit/checks")
import harness
sys.path.insert(0, "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
import camb
from camb import _ini
out = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/write-ini"
base = set(harness.parse_ini(REPO + "/fortran/eftcamb_test/parameters/base_params.ini"))
inis = sorted(glob.glob(REPO + "/fortran/eftcamb_test/parameters/[0-9]*.ini"))
tot_keys = tot_written = 0; nfail = 0
with harness.slot(1):
    for ini in inis:
        flat = harness.parse_ini(ini)
        eftkeys = [k for k in flat if k.lower().startswith("eft") or k not in base]
        p = camb.read_ini(ini)
        if p.EFTCAMB.EFTflag == 0:
            continue
        fn = os.path.join(out, os.path.basename(ini))
        camb.write_ini(p, fn, validate=False)
        w = harness.parse_ini(fn)
        wl = {k.lower() for k in w}
        written = [k for k in eftkeys if k.lower() in wl]
        tot_keys += len(eftkeys); tot_written += len(written)
        q = camb.read_ini(fn)
        if q.EFTCAMB.EFTflag != p.EFTCAMB.EFTflag: nfail += 1
        print(os.path.basename(ini), "EFTflag", p.EFTCAMB.EFTflag, "->", q.EFTCAMB.EFTflag,
              "eftkeys", len(eftkeys), "written", len(written), "dark_energy_model=", w.get("dark_energy_model"))
print("TOTAL eft-ish keys", tot_keys, "written", tot_written, "models with EFTflag lost", nfail, "of", len(inis))
