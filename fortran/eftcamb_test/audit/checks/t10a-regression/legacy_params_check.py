"""Compare the parameters recorded in the legacy (B2) *_params.ini files with the current test inis.

For each model present in both, print every key whose value differs (numeric comparison with rtol 1e-6),
restricted to keys that the current flattened ini defines.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

LEG = harness.REPO / "fortran" / "eftcamb_test" / "eftcamb_legacy" / "spectra_results"


def norm(v: str):
    v = v.strip()
    try:
        return float(v.replace("d", "e").replace("D", "e"))
    except ValueError:
        u = v.upper()
        if u in ("T", "TRUE", ".TRUE."):
            return True
        if u in ("F", "FALSE", ".FALSE."):
            return False
        return v


def main() -> None:
    skip = {"output_root", "EFTCAMB_write_background", "feedback_level", "number_of_threads", "output_file_headers"}
    for p in harness.model_inis():
        leg = LEG / f"{p.stem}_params.ini"
        if not leg.exists():
            print(f"{p.stem:34s} NO LEGACY PARAMS")
            continue
        cur = harness.parse_ini(p)
        old = {}
        for line in leg.read_text().splitlines():
            if "=" in line and not line.lstrip().startswith("#"):
                k, v = line.split("=", 1)
                old[k.strip()] = v.strip()
        diffs = []
        for k, v in cur.items():
            if k in skip or k.endswith("_file") or k.endswith("filename(1)") or k.endswith("matterpower(1)"):
                continue
            if k not in old:
                pass  # key not recorded by the legacy writer
                continue
            a, b = norm(v), norm(old[k])
            if isinstance(a, float) and isinstance(b, float):
                if abs(a - b) > 1e-6 * max(abs(a), abs(b), 1e-300):
                    diffs.append(f"{k}: cur={v} legacy={old[k]}")
            elif a != b:
                diffs.append(f"{k}: cur={v} legacy={old[k]}")
        print(f"{p.stem:34s} " + ("OK" if not diffs else "; ".join(diffs)))


if __name__ == "__main__":
    main()
