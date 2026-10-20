"""Independent check of fr-rgr-threshold-high-k-pk: skip_RGR, turn-on convergence, accuracy, GR ref, B1."""
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

OUT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
       "scratchpad/audit/work/verify/frrgr")
INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_2.ini"
CASES = {
    "A0_gr": ("A0", {"EFTflag": "0"}),
    "A0_def": ("A0", {}),
    "A0_skip": ("A0", {"EFTCAMB_skip_RGR": "T"}),
    "A0_skip_t003": ("A0", {"EFTCAMB_skip_RGR": "T", "EFTCAMB_turn_on_time": "0.003"}),
    "A0_skip_t03": ("A0", {"EFTCAMB_skip_RGR": "T", "EFTCAMB_turn_on_time": "0.03"}),
    "A0_def_acc2": ("A0", {"accuracy_boost": "2"}),
    "A0_skip_acc2": ("A0", {"EFTCAMB_skip_RGR": "T", "accuracy_boost": "2"}),
    "B1_def": ("B1", {}),
    "B1_skip": ("B1", {"EFTCAMB_skip_RGR": "T"}),
}
sel = sys.argv[1:] or list(CASES)


def go(name):
    tree, ov = CASES[name]
    return name, harness.run_model(tree, INI, f"{OUT}/{name}.npz", overrides=ov, threads=1)["status"]


with ThreadPoolExecutor(4) as ex:
    for n, s in ex.map(go, sel):
        print(n, s, flush=True)
