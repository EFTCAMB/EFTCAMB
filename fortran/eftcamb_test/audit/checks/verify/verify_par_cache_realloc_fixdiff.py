"""Check that the minimal realloc fix (two deallocate lines) leaves fresh-run outputs bit-identical to A0."""

from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

ROOT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit"
)
FIX = ROOT / "work" / "verify" / "fixmin_realloc"
OUT = ROOT / "work" / "verify" / "par_cache_realloc" / "fixdiff"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
MODELS = sys.argv[1:] or ["4_Des_F_R_LCDM_1", "5_Kmouflage_1", "2_PEFT_Omega_lin_1", "5_hdsk_TG_1"]


def one(m: str):
    out = OUT / f"{m}.npz"
    st = harness.run_model(str(FIX), str(PARAMS / f"{m}.ini"), str(out), threads=2)
    if st["status"] != "ok":
        return m, st["status"], None
    return m, "ok", harness.compare(harness.load(str(ROOT / "runs" / "A0" / f"{m}.npz")), harness.load(str(out)))


OUT.mkdir(parents=True, exist_ok=True)
with ThreadPoolExecutor(4) as ex:
    for m, status, c in ex.map(one, MODELS):
        if c is None:
            print(f"{m:24s} status={status}")
            continue
        eft = max(c["eft"].values()) if c["eft"] else 0.0
        print(f"{m:24s} bit_identical={c['bit_identical']} worst_cl={c['worst_cl']} worst_pk={c['worst_pk']} "
              f"sigma8={c['sigma8']} max_eft={eft:.2e}")
