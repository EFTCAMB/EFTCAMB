"""Verifier (read-parameters-roundtrip-bytes): independent check of camb/eftcamb.py EFTCAMB.read_parameters.

Checks, each in a fresh subprocess (ERROR STOP in the Fortran ini reader kills the interpreter):
  1. raw: call the Fortran getter (get_read_parameters) directly via ctypes and decode the buffers ourselves, to see
     what Fortran stores (attribution: Fortran side vs python autoconvert).
  2. pydict: python-only path (no read_ini): initialize_parameters with a dict where EFTB0 is the *string* '1.d-4'
     (and output_root='abc'); Fortran accepts it (param_values -> 1e-4) but read_parameters returns bytes.
  3. census: over all test inis, camb.read_ini + read_parameters: which keys come back as bytes, and which of those are
     keys that EFTCAMB itself reads (keys present after a python re-initialization).
  4. roundtrip (as shipped): initialize_parameters(read_parameters()) per ini; record crash / garbled keys.
  5. roundtrip (decoded): same with a locally patched read_parameters that decodes and maps Fortran D exponents
     (the proposed fix, monkeypatched in this process only), comparing EFT ctypes fields, param_values and
     read parameters before/after.

usage: python readparams_roundtrip_indep.py [tree] [outjson]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

REPO = HERE.parents[4]
FORT = REPO / "fortran"
OUT_DEFAULT = Path(
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/work/verify/readparams_rt/indep.json"
)

CHILD_COMMON = r"""
import sys, json, re, copy
sys.path.insert(0, sys.argv[1])
import camb
from camb import eftcamb as E
from ctypes import c_int, byref

def raw_read(eft):
    n = c_int(0)
    k = (' ' * 200).encode(); v = (' ' * 200).encode()
    E.EFTCAMB_get_read_parameters(byref(eft), k, v, byref(n), byref(c_int(0)))
    out = {}
    for i in range(n.value):
        k = (' ' * 200).encode(); v = (' ' * 200).encode()
        E.EFTCAMB_get_read_parameters(byref(eft), k, v, byref(c_int(i + 1)), byref(c_int(1)))
        out[k.decode().strip()] = v.decode().strip()
    return out

DEXP = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)[dD][+-]?\d+$')
def conv(s):
    if s == 'T':
        return True
    if s == 'F':
        return False
    for fn in (int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    if DEXP.match(s):
        return float(s.replace('d', 'e').replace('D', 'E'))
    return s

def fixed_read_parameters(self):
    return {k: conv(v) for k, v in raw_read(self).items()}

def fields(eft):
    return {f[0]: getattr(eft, f[0]) for f in E.EFTCAMB._fields_}
"""

CHILD_RAW = CHILD_COMMON + r"""
p = camb.read_ini(sys.argv[2])
raw = raw_read(p.EFTCAMB)
d = p.EFTCAMB.read_parameters()
print(json.dumps({'raw_EFTB0': raw.get('EFTB0'), 'raw_output_root': raw.get('output_root'),
                  'rp_EFTB0': repr(d.get('EFTB0')), 'rp_output_root': repr(d.get('output_root')),
                  'param_values': p.EFTCAMB.param_values(), 'param_names': p.EFTCAMB.param_names()}))
"""

CHILD_PYDICT = CHILD_COMMON + r"""
q = camb.CAMBparams(); q.set_cosmology(H0=67)
q.EFTCAMB.initialize_parameters(q, {'EFTflag': 3, 'DesignerEFTmodel': 1, 'EFTwDE': 0, 'EFTB0': sys.argv[2],
                                    'output_root': 'abc'}, print_header=False)
d = q.EFTCAMB.read_parameters()
from camb.eftcamb_tree import check_params
name, inp, outp = check_params({'EFTflag': 3, 'DesignerEFTmodel': 1, 'EFTwDE': 0, 'EFTB0': 0.5, 'output_root': 'abc'})
print(json.dumps({'EFTB0_input': sys.argv[2], 'rp_EFTB0': repr(d.get('EFTB0')), 'rp_output_root': repr(d.get('output_root')),
                  'param_values': q.EFTCAMB.param_values(), 'keys_read': sorted(d),
                  'check_params_output_root_in_input_dict': 'output_root' in inp,
                  'check_params_output_root_in_output_dict': repr(outp.get('output_root'))}))
"""

CHILD_CENSUS = CHILD_COMMON + r"""
res = {}
for ini in sys.argv[2:]:
    p = camb.read_ini(ini)
    d = p.EFTCAMB.read_parameters()
    q = camb.CAMBparams(); q.set_cosmology(H0=67)
    q.EFTCAMB.initialize_parameters(q, fixed_read_parameters(p.EFTCAMB), print_header=False)
    eft_keys = set(raw_read(q.EFTCAMB))
    bk = sorted(k for k, v in d.items() if isinstance(v, bytes))
    res[ini.split('/')[-1]] = {'n_keys': len(d), 'bytes_keys': bk, 'bytes_eft_keys': sorted(set(bk) & eft_keys),
                               'dexp_eft_keys': sorted(k for k in set(bk) & eft_keys
                                                       if DEXP.match(d[k].decode()))}
print(json.dumps(res))
"""

CHILD_RT = CHILD_COMMON + r"""
mode = sys.argv[3]
p = camb.read_ini(sys.argv[2])
if mode == 'fixed':
    E.EFTCAMB.read_parameters = fixed_read_parameters
d1 = p.EFTCAMB.read_parameters()
f1, v1 = fields(p.EFTCAMB), p.EFTCAMB.param_values()
q = camb.CAMBparams(); q.set_cosmology(H0=67)
q.EFTCAMB.initialize_parameters(q, d1, print_header=False)
d2 = q.EFTCAMB.read_parameters()
f2, v2 = fields(q.EFTCAMB), q.EFTCAMB.param_values()
diff_keys = sorted(k for k in d2 if d2[k] != d1.get(k))
print('RTJSON' + json.dumps({'fields_equal': f1 == f2, 'values_equal': v1 == v2, 'n_values': len(v1),
                             'changed_read_keys': {k: [repr(d1.get(k)), repr(d2[k])] for k in diff_keys}}))
"""


def run(code: str, args: list[str], tree: str) -> subprocess.CompletedProcess:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    return subprocess.run([sys.executable, "-c", code, tree, *args], env=env, capture_output=True, text=True,
                          cwd=str(FORT), timeout=600)


def last_json(r: subprocess.CompletedProcess):
    for ln in reversed(r.stdout.splitlines()):
        if ln.startswith("RTJSON"):
            return json.loads(ln[6:])
        if ln.startswith("{"):
            return json.loads(ln)
    return None


if __name__ == "__main__":
    tree = str(harness.tree_path(sys.argv[1] if len(sys.argv) > 1 else "A0"))
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else OUT_DEFAULT
    out.parent.mkdir(parents=True, exist_ok=True)
    inis = [str(p) for p in harness.model_inis()]
    report = {"tree": tree, "n_inis": len(inis)}
    with harness.slot(1):
        r = run(CHILD_RAW, [str(harness.PARAMS_DIR / "4_Des_F_R_LCDM_2.ini")], tree)
        report["raw_4_Des_F_R_LCDM_2"] = last_json(r) or {"rc": r.returncode, "tail": r.stderr[-500:]}
        print("raw:", report["raw_4_Des_F_R_LCDM_2"])
        report["pydict"] = []
        for s in ("1.d-4", "1.0D-04", "1e-4"):
            r = run(CHILD_PYDICT, [s], tree)
            report["pydict"].append(last_json(r) or {"rc": r.returncode, "tail": r.stderr[-500:]})
            print("pydict:", report["pydict"][-1])
        r = run(CHILD_CENSUS, inis, tree)
        census = last_json(r) or {}
        report["census_rc"] = r.returncode
        report["census"] = census
        n_bytes_root = sum("output_root" in c["bytes_eft_keys"] for c in census.values())
        bytes_eft = sorted({k for c in census.values() for k in c["bytes_eft_keys"]})
        dexp = {m: c["dexp_eft_keys"] for m, c in census.items() if c["dexp_eft_keys"]}
        report["census_summary"] = {"inis": len(census), "output_root_bytes_eft": n_bytes_root,
                                    "all_bytes_eft_keys": bytes_eft, "dexp_eft_keys": dexp}
        print("census:", report["census_summary"])
        for mode in ("shipped", "fixed"):
            rows = {}
            for ini in inis:
                r = run(CHILD_RT, [ini, mode], tree)
                j = last_json(r)
                rows[Path(ini).stem] = {"rc": r.returncode, **(j or {"tail": (r.stdout + r.stderr)[-300:]})}
            ok = [m for m, v in rows.items() if v["rc"] == 0 and v.get("fields_equal") and v.get("values_equal")
                  and not v.get("changed_read_keys")]
            crash = {m: v.get("tail", "").strip().splitlines()[-2:] for m, v in rows.items() if v["rc"] != 0}
            garbled = {m: v["changed_read_keys"] for m, v in rows.items() if v["rc"] == 0 and v.get("changed_read_keys")}
            report[f"roundtrip_{mode}"] = {"identical": len(ok), "crash": crash,
                                           "garbled_count": len(garbled),
                                           "garbled_example": dict(list(garbled.items())[:2]),
                                           "rows": rows}
            print(f"roundtrip {mode}: identical {len(ok)}/{len(inis)}, crash {len(crash)} {crash}, "
                  f"garbled {len(garbled)}; example {dict(list(garbled.items())[:1])}")
    out.write_text(json.dumps(report, indent=1, default=str))
    print("wrote", out)
