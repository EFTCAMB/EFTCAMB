"""(vi) camb.set_params positivity-bounds allowlist: model names vs allowlist, and behaviour of the guard."""
import json, os, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import harness

CODE = r'''
import sys, json
sys.path.insert(0, TREE)
import camb
out = {}
base = dict(H0=70, ombh2=0.0226, omch2=0.112, dark_energy_model="EFTCAMB", feedback_level=0)
models = {
  "PureEFT": dict(EFTflag=1, PureEFTmodel=1, EFTwDE=0, PureEFTmodelGamma1=1, EFTGamma10=0.1),
  "RPH": dict(EFTflag=2, AltParEFTmodel=1, RPHwDE=0, RPHmassPmodel=1, RPHkineticitymodel=1, RPHbraidingmodel=0,
              RPHtensormodel=0, RPHmassP0=0.1, RPHkineticity0=0.1),
  "OLgamma": dict(EFTflag=2, AltParEFTmodel=2),
  "ShiftSym": dict(EFTflag=2, AltParEFTmodel=4),
  "fR": dict(EFTflag=3, DesignerEFTmodel=1, EFTwDE=0, EFTB0=0.01),
  "3p1": dict(EFTflag=3, DesignerEFTmodel=3),
  "Kmouflage": dict(EFTflag=4, FullMappingEFTmodel=3),
  "Quintessence": dict(EFTflag=4, FullMappingEFTmodel=4),
  "SCG": dict(EFTflag=4, FullMappingEFTmodel=6),
}
for name, m in models.items():
    if name != ONLY:
        continue
    rec = {}
    for pb in (None, True, False):
        kw = dict(base, **m)
        if pb is not None:
            kw["EFT_positivity_bounds"] = pb
        try:
            cp = camb.set_params(**kw)
            rec[str(pb)] = "ok"
            rec["model_name"] = cp.EFTCAMB.model_name()
            rec["field_EFT_positivity_bounds_" + str(pb)] = bool(cp.EFTCAMB.EFT_positivity_bounds)
        except Exception as e:
            rec[str(pb)] = type(e).__name__ + ": " + str(e)[:160]
    out[name] = rec
print("T6JSON " + json.dumps(out))
'''
tree = str(harness.tree_path(sys.argv[1] if len(sys.argv) > 1 else "A0"))
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
for only in ("PureEFT", "RPH", "OLgamma", "ShiftSym", "fR", "3p1", "Kmouflage", "Quintessence", "SCG"):
    with harness.slot(1):
     p = subprocess.run([sys.executable, "-c", "TREE=%r\nONLY=%r\n" % (tree, only) + CODE], capture_output=True,
                       text=True, env=env, cwd=str(Path(__file__).resolve().parent), timeout=600)
    got = [json.loads(l[7:]) for l in p.stdout.splitlines() if l.startswith("T6JSON ")]
    print(only, "rc", p.returncode, json.dumps(got[0] if got else None),
          "" if got else (p.stdout[-300:] + p.stderr[-300:]).replace("\n", " | "))
