"""Track 9 (ii)+(iii): read_parameters round trip and set_params equivalence for every test ini (tree A0).

Per model, in one subprocess:
  R1 = read_ini(model.ini)                                    (Fortran ini path, reference)
  d1 = R1.EFTCAMB.read_parameters()
  R2 = read_ini(GR wrapper: DEFAULT(model.ini) + EFTflag = 0) then EFTCAMB.initialize_parameters(R2, d1')   (ii)
  R3 = read_ini(GR wrapper) then camb.set_params(cp=R3, **d1')                                               (iii)
  R0 = read_ini(model.ini) again, only if R2 or R3 is not bit identical (determinism control)
where d1' = d1 with bytes values decoded to str (raw bytes values make Fortran 'error stop', see repro).
Also checks that read_parameters after the round trip returns the same dict (values and Python types).

Usage: python roundtrip.py [glob] [--jobs N]
"""

from __future__ import annotations

import argparse
import fnmatch
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from t9lib import WORK, model_inis, run_snippet

OUT = WORK / "roundtrip"

CODE = r"""
import time
model = PARAMS / {model!r}
out = WORK / 'roundtrip'
FAST = 'l_max_scalar = 2200\nl_max_tensor = 1500\nk_eta_max_tensor = 3000\n'
fast = out / ({stem!r} + '.fast.ini')
fast.write_text(f'DEFAULT({{model}})\n' + FAST)
gr = out / ({stem!r} + '.gr.ini')
gr.write_text(f'DEFAULT({{model}})\nEFTflag = 0\n' + FAST)
model = fast

def summary(pars):
    r = camb.get_results(pars)
    cls = r.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)
    res = {{'total': np.asarray(cls['total']), 'unlensed_scalar': np.asarray(cls['unlensed_scalar'])}}
    if pars.WantTransfer:
        k, z, pk = r.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
        res['pk'] = np.asarray(pk); res['sigma8'] = np.atleast_1d(r.get_sigma8())
    return res

def same(a, b):
    ident = all(np.array_equal(a[k], b[k]) for k in a)
    rel = 0.0
    for k in a:
        x, y = a[k], b[k]
        if k in ('total', 'unlensed_scalar'):
            x, y = x[2:2500, 0], y[2:2500, 0]
        m = y != 0
        if m.any():
            rel = max(rel, float(np.max(np.abs(x[m] / y[m] - 1))))
    return ident, rel

def clean(d):
    return {{k: (v.decode() if isinstance(v, bytes) else v) for k, v in d.items()}}

t0 = time.time()
p1 = camb.read_ini(str(model))
d1 = p1.EFTCAMB.read_parameters()
emit(EFTflag=int(p1.EFTCAMB.EFTflag), model_name=p1.EFTCAMB.model_name(), nkeys=len(d1),
     types={{k: type(v).__name__ for k, v in d1.items()}},
     bytes_values={{k: repr(v) for k, v in d1.items() if isinstance(v, bytes)}})
R1 = summary(p1)
d1c = clean(d1)

# (ii) dictionary round trip on a GR copy of the same ini
p2 = camb.read_ini(str(gr))
assert p2.EFTCAMB.EFTflag == 0
p2.EFTCAMB.initialize_parameters(p2, d1c, print_header=False)
d2 = p2.EFTCAMB.read_parameters()
d2c = clean(d2)
dict_same = (d2c == d1c) and all(type(d2c[k]) is type(d1c[k]) for k in d1c)
diff = {{k: [repr(d1c.get(k)), repr(d2c.get(k))] for k in set(d1c) | set(d2c) if d1c.get(k) != d2c.get(k)
        or type(d1c.get(k)) is not type(d2c.get(k))}}
R2 = summary(p2)
i2, rel2 = same(R2, R1)
emit(dict_roundtrip_equal=dict_same, dict_diff=diff, ii_identical=i2, ii_maxrel=rel2)

# (iii) camb.set_params(cp=..., **eft_keys)
p3 = camb.read_ini(str(gr))
camb.set_params(cp=p3, **d1c)
d3 = p3.EFTCAMB.read_parameters()
R3 = summary(p3)
i3, rel3 = same(R3, R1)
emit(iii_identical=i3, iii_maxrel=rel3, iii_dict_equal=(clean(d3) == d1c))
if not (i2 and i3):
    p0 = camb.read_ini(str(model))
    R0 = summary(p0)
    i0, rel0 = same(R0, R1)
    emit(control_identical=i0, control_maxrel=rel0)
emit(elapsed=time.time() - t0)
"""


def one(ini: Path) -> dict:
    r = run_snippet(CODE.format(model=ini.name, stem=ini.stem), timeout=1500)
    rec = {"model": ini.stem, "status": r["status"], **r["result"]}
    if r["status"] != "ok":
        rec["tail"] = r["tail"][-10:]
    (OUT / f"{ini.stem}.json").write_text(json.dumps(rec, indent=1))
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern", nargs="?", default="*")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    inis = [p for p in model_inis() if fnmatch.fnmatch(p.stem, a.pattern)]
    if a.skip_done:
        inis = [p for p in inis if not (OUT / f"{p.stem}.json").exists()]
    with ThreadPoolExecutor(a.jobs) as ex:
        for rec in ex.map(one, inis):
            print(f"{rec['model']:32s} {rec['status']:10s} dict={rec.get('dict_roundtrip_equal')} "
                  f"ii={rec.get('ii_identical')} ({rec.get('ii_maxrel', float('nan')):.1e}) "
                  f"iii={rec.get('iii_identical')} ({rec.get('iii_maxrel', float('nan')):.1e}) "
                  f"ctrl={rec.get('control_identical', '-')} bytes={list(rec.get('bytes_values', {}))}", flush=True)


if __name__ == "__main__":
    main()
