"""Verifier (t3-eft-turnon-before-tauminn-lmax): scan l_max_scalar (and transfer-only runs) for 5_ADE_1 on a tree,
collect pass/fail and the tauminn printed in the EFTCAMB refusal message, to see how tauminn scales.

Usage: python tauminn_lmax_scan.py <tree> <tag> lmax1 lmax2 ... [--transfer-only]
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "5_ADE_1.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "tauminn"


def job(tree: str, tag: str, case: str, overrides: dict) -> dict:
    out = OUT / f"{tag}_{case}.npz"
    st = harness.run_model(tree, INI, out, overrides=overrides, threads=1, timeout=900)
    log = out.with_suffix(".log").read_text(errors="replace")
    msg = re.findall(r"EFTCAMB starts before thermo tauminn[^\n]*= *([^,]+), *([^,]+), *([^\s,]+?)\.?\s*$", log, re.M)
    return {"case": case, "status": st["status"], "msg": msg[:1]}


def main() -> None:
    tree, tag = sys.argv[1], sys.argv[2]
    cases = {}
    for a in sys.argv[3:]:
        if a == "--transfer-only":
            cases["transfer_only"] = {"get_scalar_cls": "F", "get_tensor_cls": "F", "get_vector_cls": "F",
                                      "do_lensing": "F"}
        elif a.startswith("k"):
            cases[f"transfer_only_kmax{a[1:]}"] = {"get_scalar_cls": "F", "get_tensor_cls": "F",
                                                   "do_lensing": "F", "transfer_kmax": a[1:]}
        else:
            cases[f"lmax{a}"] = {"l_max_scalar": a, "get_tensor_cls": "F"}
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(lambda kv: job(tree, tag, kv[0], kv[1]), cases.items()))
    for r in res:
        print(json.dumps(r))


if __name__ == "__main__":
    main()
