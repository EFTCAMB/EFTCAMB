"""Code verdict (get_background with stability on) for all test models on a given tree."""
import sys, json
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, ".")
import t6lib, harness
tree = sys.argv[1]
extra = dict(s.split("=", 1) for s in sys.argv[2:])
inis = harness.model_inis()
def one(ini):
    r = t6lib.probe(tree, ini, overrides={"feedback_level": "0", **extra}, name=ini.stem)
    return ini.stem, r["status"], r.get("bg_error")
with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(one, inis))
bad = [x for x in res if x[1] != "ok" or x[2] is not None]
print(json.dumps({"tree": tree, "n": len(res), "n_unstable_or_failed": len(bad), "bad": bad}, indent=1))
