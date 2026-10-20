"""Code verdict (get_background, stability forced on) for all test models on a tree, per-model timeout, all results."""
import sys, json
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, ".")
import t6lib, harness
tree = sys.argv[1]
extra = dict(s.split("=", 1) for s in sys.argv[2:])
def one(ini):
    r = t6lib.probe(tree, ini, overrides={"feedback_level": "2", **extra}, name=ini.stem, timeout=150)
    nanmsg = sorted({l.strip() for l in r["log"] if " is NaN" in l})[:3]
    first = next((l.strip() for l in r["log"] if " is NaN" in l), None)
    return {"model": ini.stem, "status": r["status"],
            "verdict": "stable" if r.get("bg_error") is None and r["status"] == "ok" else "UNSTABLE/FAIL",
            "err": (r.get("bg_error") or "")[-45:], "first_nan": first}
with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(one, harness.model_inis()))
print(json.dumps(res, indent=0))
