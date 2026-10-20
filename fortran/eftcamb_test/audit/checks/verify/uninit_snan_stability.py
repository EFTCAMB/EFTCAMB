"""Run models on flags-only snan tree with EFTCAMB_feedback_level=2; report status and first NaN field."""
import re, sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
W = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/uninit"
tree = sys.argv[1] if len(sys.argv) > 1 else f"{W}/snan"
tag = sys.argv[2] if len(sys.argv) > 2 else "snan"
models = ["4_Des_F_R_LCDM_1", "5_quint_exponential", "5_Kmouflage_3", "5_ADE_1", "5_Horava_eta_1", "5_Kmimic_1",
          "5_scg", "5_hdsk_TG_1", "1_Pure_EFT_std", "2_Lambda_alpha"]
cheap = {"get_transfer": "F", "get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "100",
         "EFTCAMB_feedback_level": "2"}
for m in models:
    ini = harness.PARAMS_DIR / f"{m}.ini"
    if not ini.exists():
        cands = sorted(harness.PARAMS_DIR.glob(m.split("_")[0] + "_*" + m.split("_", 1)[1] + "*.ini"))
        if not cands:
            print(m, "no ini"); continue
        ini = cands[0]
    st = harness.run_model(tree, ini, f"{W}/{tag}_{ini.stem}.npz", overrides=cheap, threads=1, timeout=500)
    log = open(st["log"]).read()
    nans = re.findall(r"^\s*(\S+)\s+is NaN", log, re.M)
    unstable = "unstable" in log.lower()
    print(f"{ini.stem:24s} {st['status']:12s} unstable_msg={unstable} first_NaN={nans[:1]} n_NaN_lines={len(nans)}")
