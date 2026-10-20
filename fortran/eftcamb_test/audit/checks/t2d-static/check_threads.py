"""Thread-race probe: run selected models on A0 with OMP_NUM_THREADS=1 and =6 and compare (bit identity expected,
since every k mode is independent and the k sums are done in a fixed order).

Usage: python check_threads.py OUTDIR [model ...]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

MODELS = ["2_PEFT_gamma3_const_1", "2_PEFT_Omega_power_1", "3_RPH_KT_power_1", "4_Des_F_R_wCDM_1",
          "4_Des_mc5e_taylor_1", "5_hdsk_TG_2", "5_quint_axion", "5_Kmouflage_1", "5_ADE_1", "5_Horava_eta_1",
          "5_scg", "2_PEFT_wDE_wCDM_1"]


def main():
    out = Path(sys.argv[1])
    models = sys.argv[2:] or MODELS
    res = {}
    for m in models:
        ini = harness.PARAMS_DIR / f"{m}.ini"
        st1 = harness.run_model("A0", ini, out / f"{m}_t1.npz", threads=1)
        st6 = harness.run_model("A0", ini, out / f"{m}_t6.npz", threads=6)
        if st1["status"] != "ok" or st6["status"] != "ok":
            res[m] = {"status": [st1["status"], st6["status"]]}
        else:
            c = harness.compare(harness.load(out / f"{m}_t1.npz"), harness.load(out / f"{m}_t6.npz"))
            res[m] = {"bit_identical": c["bit_identical"], "worst_cl": c["worst_cl"], "worst_pk": c["worst_pk"],
                      "sigma8": c["sigma8"], "eft": max(c["eft"].values()) if c["eft"] else 0.0}
        print(m, json.dumps(res[m]), flush=True)
    (out / "threads_summary.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
