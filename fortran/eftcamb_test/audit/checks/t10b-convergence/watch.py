"""Print one line per finished run (json status file) under the given run dirs; exits when n_expected are done."""
import json
import sys
import time
from pathlib import Path

dirs = [Path(d) for d in sys.argv[2:]]
n_expected = int(sys.argv[1])
seen = set()
while True:
    for d in dirs:
        for f in sorted(d.glob("*.json")):
            if f not in seen:
                seen.add(f)
                try:
                    j = json.loads(f.read_text())
                    print(f"{d.name}/{f.stem} {j['status']} {j['elapsed']}s", flush=True)
                except Exception as e:  # noqa: BLE001
                    print(f"{d.name}/{f.stem} unreadable {e}", flush=True)
    if len(seen) >= n_expected:
        break
    time.sleep(20)
