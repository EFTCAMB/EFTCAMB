"""Independent check of finding write-ini-drops-eftcamb (verify track).

Differences from the finder's repro:
  * uses the module-level camb.write_ini (not CAMBparams.write_ini) for the unvalidated write;
  * the EFT key set is taken from the ini TEXT (flattened model ini minus base_params keys, plus every EFT* key),
    not from EFTCAMB.read_parameters();
  * written files are parsed with an independent minimal parser;
  * root-cause test: monkeypatch camb._ini._update_ini_state_from_params to append the EFT keys (values copied
    verbatim from the original ini text) and rerun write_ini(validate=True) -> does validation then pass and do
    spectra round-trip?
Each model runs in its own child process (sequentially) inside one harness CPU slot.

Usage: python writeini_indep.py [--tree PATH] [--models GLOB] [--spectra N] --out JSON
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

CHILD = r"""
import sys, os, json
cfg = json.loads(sys.argv[1])
sys.path.insert(0, cfg['tree'])
import numpy as np
import camb
from camb import _ini
from camb.baseconfig import CAMBValueError

res = {'model': cfg['name']}
tmp = cfg['tmp']
fast = os.path.join(tmp, 'model.ini')
open(fast, 'w').write(f"DEFAULT({cfg['ini']})\nl_max_scalar = 1200\nget_tensor_cls = F\nget_transfer = F\n"
                      f"EFTCAMB_write_background = F\noutput_root = {tmp}/out\n")
camb.config.print_fortran_warnings = False
p = camb.read_ini(fast)
res['EFTflag'] = int(p.EFTCAMB.EFTflag)
res['camb_file'] = camb.__file__

def keys_of(path):
    out = {}
    for line in open(path):
        s = line.split('#', 1)[0].strip()
        if '=' in s:
            k, v = s.split('=', 1)
            out[k.strip()] = v.strip()
    return out

# 1) unvalidated write through the module-level function
w0 = os.path.join(tmp, 'written_novalidate.ini')
camb.write_ini(p, w0, validate=False)
wk = keys_of(w0)
eft_keys = cfg['eft_keys']
res['n_eft_keys_text'] = len(eft_keys)
res['n_eft_keys_written'] = sum(k in wk for k in eft_keys)
res['eft_keys_written'] = [k for k in eft_keys if k in wk]
res['EFTflag_in_written'] = 'EFTflag' in wk
q = camb.read_ini(w0)
res['reread_EFTflag'] = int(q.EFTCAMB.EFTflag)

# 2) default validated write
try:
    p.write_ini(os.path.join(tmp, 'written_validate.ini'))
    res['validate_error'] = None
except CAMBValueError as e:
    res['validate_error'] = str(e).split(': ', 1)[-1][:200]

# 3) root-cause test: append EFT keys (verbatim ini values, only the ones EFTCAMB actually read)
read_names = set(p.EFTCAMB.read_parameters().keys())
flat = cfg['flat']
add = {k: v for k, v in flat.items() if k in read_names and k != 'output_root'}
res['n_added'] = len(add)
orig = _ini._update_ini_state_from_params
def patched(params, state):
    orig(params, state)
    for k, v in add.items():
        state.set(k, v)
_ini._update_ini_state_from_params = patched
w1 = os.path.join(tmp, 'written_patched.ini')
try:
    p.write_ini(w1)
    res['patched_validate_error'] = None
except CAMBValueError as e:
    res['patched_validate_error'] = str(e).split(': ', 1)[-1][:200]
r = camb.read_ini(w1)
res['patched_reread_EFTflag'] = int(r.EFTCAMB.EFTflag)
rp0 = {k: str(v) for k, v in p.EFTCAMB.read_parameters().items() if k != 'output_root'}
rp1 = {k: str(v) for k, v in r.EFTCAMB.read_parameters().items() if k != 'output_root'}
res['patched_read_parameters_equal'] = rp0 == rp1
if rp0 != rp1:
    res['patched_read_parameters_diff'] = {k: (rp0.get(k), rp1.get(k)) for k in set(rp0) | set(rp1)
                                           if rp0.get(k) != rp1.get(k)}

if cfg['spectra']:
    def tt(pp):
        return camb.get_results(pp).get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:1000, 0]
    c0, cq, cr = tt(p), tt(q), tt(r)
    res['TT_reread_unpatched_maxrel'] = float(np.max(np.abs(cq / c0 - 1)))
    res['TT_reread_patched_maxrel'] = float(np.max(np.abs(cr / c0 - 1)))
print('RESULT_JSON ' + json.dumps(res))
"""


def eft_keys_from_text(ini: Path) -> tuple[list[str], dict[str, str]]:
    flat = harness.parse_ini(ini)
    base = harness.parse_ini(harness.PARAMS_DIR / "base_params.ini")
    hdsk = harness.parse_ini(harness.PARAMS_DIR / "hdsk_base_params.ini")
    camb_keys = {k for k in list(base) + list(hdsk) if not k.startswith("EFT")}
    keys = sorted(k for k in flat if k != "output_root" and (k.startswith("EFT") or k not in camb_keys))
    return keys, flat


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=str(harness.TREES["A0"]))
    ap.add_argument("--models", default="*.ini")
    ap.add_argument("--spectra", default="", help="comma list of model names for which to compute TT")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    spectra = set(filter(None, a.spectra.split(",")))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    results = []
    with harness.slot(1):
        for ini in harness.model_inis(a.models):
            keys, flat = eft_keys_from_text(ini)
            # same overrides as the child's wrapper ini
            flat.update({"l_max_scalar": "1200", "get_tensor_cls": "F", "get_transfer": "F",
                         "EFTCAMB_write_background": "F"})
            with tempfile.TemporaryDirectory() as tmp:
                cfg = dict(tree=a.tree, ini=str(ini), name=ini.stem, tmp=tmp, eft_keys=keys, flat=flat,
                           spectra=ini.stem in spectra)
                r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, json.dumps(cfg)], cwd=tmp,
                                   env=env, capture_output=True, text=True, timeout=900)
            line = [x for x in r.stdout.splitlines() if x.startswith("RESULT_JSON ")]
            if line:
                res = json.loads(line[-1][len("RESULT_JSON "):])
            else:
                res = {"model": ini.stem, "child_error": f"rc={r.returncode}", "stderr": r.stderr[-600:]}
            results.append(res)
            print(json.dumps({k: v for k, v in res.items() if k not in ("eft_keys_written", "camb_file")}),
                  flush=True)
    Path(a.out).write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
