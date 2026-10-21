"""Run K-mouflage / K-mimic with EFT_positivity_bounds=T on A0 to see whether the positivity check is reached."""
import sys
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
T = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0"
O = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/hdddot_pos"
inis = {m: i for m, i in ((p.stem, p) for p in harness.model_inis())}
for m in ("5_Kmouflage_1", "5_Kmimic_1"):
    st = harness.run_model(T, str(inis[m]), f"{O}/{m}.npz", overrides={"EFT_positivity_bounds": "T"}, threads=1, timeout=500)
    print(m, st["status"], st.get("meta_errors") if "meta_errors" in st else "")
    print("\n".join(l for l in st["log_tail"].splitlines()[-8:]))
