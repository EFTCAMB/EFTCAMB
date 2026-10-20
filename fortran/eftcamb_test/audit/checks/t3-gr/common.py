"""Shared helpers for track t3-gr (GR / LCDM limits)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

AUDIT = harness.AUDIT_ROOT
WORK = AUDIT / "work" / "t3-gr"
RUNS = WORK / "runs"
PARAMS = harness.PARAMS_DIR
GR_INI = PARAMS / "1_EFT_GR.ini"
INIS = HERE.parents[1] / "inis"


def summary(res: dict) -> dict:
    """Compact summary of a harness.compare result."""
    out = {"bit": res["bit_identical"], "worst_cl": res["worst_cl"], "worst_pk": res["worst_pk"]}
    for key, d in res["cl"].items():
        for s, v in d.items():
            out[f"{key}.{s}"] = v
    for key in ("sigma8", "bg_H", "bg_tau"):
        if key in res:
            out[key] = res[key]
    return out


def fmt(res: dict) -> str:
    s = summary(res)
    tot = res["cl"].get("total", {})
    return (f"bit={s['bit']!s:5} worst_cl={s['worst_cl']:.2e} worst_pk={s['worst_pk']:.2e} "
            + " ".join(f"{k}={v:.1e}" for k, v in tot.items())
            + f" s8={s.get('sigma8', float('nan')):.1e} H={s.get('bg_H', float('nan')):.1e}")


def run(tree, name, ini=GR_INI, overrides=None, python_set=None, sub="", **kw) -> Path:
    out = RUNS / sub / f"{name}.npz" if sub else RUNS / f"{name}.npz"
    st = harness.run_model(tree, ini, out, overrides=overrides, python_set=python_set, **kw)
    if st["status"] != "ok":
        print(name, st["status"], "\n".join(st.get("log_tail", [])[-15:]))
    return out


def dump(obj, path) -> None:
    Path(path).write_text(json.dumps(obj, indent=2, default=float))
