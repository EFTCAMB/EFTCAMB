"""Independent attribution of the EFT_IC_type 2/3 finding using an instrumented copy of A0 (env-var switches).

IC23_PI_MODE=attractor replaces only the pi, pi' initial values of IC 2/3 by the IC-1 (quasi-static attractor) formula,
IC23_ADIAB=gr replaces only the CAMB adiabatic ICs by the GR ones; IC23_DUMP prints series vs attractor at turn-on.

usage: python ic23_attrib.py <model> <tag> [key=value ...]   (runs the standard variant set for that model)
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work" / "verify" / "ic23"
SW = W / "A0_sw"
PAR = harness.PARAMS_DIR

model, tag = sys.argv[1], sys.argv[2]
extra = dict(a.split("=", 1) for a in sys.argv[3:])
ic = extra.pop("IC", "2")
variants = {
    "ic1": ({"EFT_IC_type": "1"}, {}),
    f"ic{ic}_series": ({"EFT_IC_type": ic}, {"IC23_DUMP": "1"}),
    f"ic{ic}_attr": ({"EFT_IC_type": ic}, {"IC23_PI_MODE": "attractor"}),
    f"ic{ic}_grad": ({"EFT_IC_type": ic}, {"IC23_ADIAB": "gr"}),
    f"ic{ic}_attr_grad": ({"EFT_IC_type": ic}, {"IC23_PI_MODE": "attractor", "IC23_ADIAB": "gr"}),
}
out = W / "attrib" / f"{model}_{tag}"
out.mkdir(parents=True, exist_ok=True)


def one(item):
    name, (ov, env) = item
    ov = dict(ov, **extra)
    return name, harness.run_model(SW, PAR / f"{model}.ini", out / f"{name}.npz", overrides=ov, extra_env=env,
                                   threads=1, timeout=1500)


with ThreadPoolExecutor(4) as ex:
    st = dict(ex.map(one, variants.items()))

res = {}
ref = harness.load(out / "ic1.npz") if st["ic1"]["status"] == "ok" else None
for name in variants:
    s = st[name]
    row = {"status": s["status"]}
    if s["status"] == "ok" and ref is not None:
        d = harness.load(out / f"{name}.npz")
        row["errors"] = d["meta"].get("errors")
        try:
            c = harness.compare(d, ref)
            cl = c["cl"]["total"]
            row.update(TT=cl["TT"], EE=cl["EE"], TE=cl["TE"], pk0=c["pk"][0] if c["pk"] else None,
                       sigma8=[float(x) for x in d["sigma8"]], sigma8_ref=[float(x) for x in ref["sigma8"]])
        except Exception as e:  # noqa: BLE001
            row["cmp_error"] = repr(e)
    else:
        row["log_tail"] = s.get("log_tail", [])[-4:]
    res[name] = row
(out / "summary.json").write_text(json.dumps(res, indent=1, default=float))
print(json.dumps(res, indent=1, default=float))
