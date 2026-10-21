"""Track 9 (vi): F-001 follow-up. Enumerate the character-argument interfaces called from python, and test the
proposed fix (pass gfortran's hidden string length, size_t, at the end of the argument list, into a mutable
create_string_buffer) on the -fcheck=all tree AD, where the current calls abort.

Each step runs in its own subprocess:
  current_<getter>  : the wrapper method as shipped (expect a -fcheck abort on AD, fine on A0)
  fixed_<getter>    : direct ctypes call with explicit hidden lengths (expect correct values on AD and A0)
Also checks (A0) that the shipped python passes immutable bytes objects that Fortran overwrites, and the
multibyte-name decode failure.

Usage: python probe_f001.py [tree]   (default AD)
"""

from __future__ import annotations

import json
import sys

from t9lib import run_snippet

MODEL = "2_PEFT_Omega_power_1.ini"
SETUP = f"""
from ctypes import byref, c_int, c_size_t, create_string_buffer
from camb import eftcamb as E
p = camb.read_ini(str(PARAMS / {MODEL!r}), no_validate=True)
eft = p.EFTCAMB
"""

STEPS = {
    "current_model_name": "emit(value=eft.model_name())",
    "current_param_names": "emit(value=eft.param_names())",
    "current_param_labels": "emit(value=eft.param_labels())",
    "current_read_parameters": "d = eft.read_parameters(); emit(n=len(d), first=list(d.items())[:3])",
    "current_set_model_name": "eft.set_model_name('renamed'); emit(done=True)",
    "fixed_model_name": """
buf = create_string_buffer(200)
E.EFTCAMB_get_model_name(byref(eft), buf, c_size_t(200))
emit(value=buf.raw.rstrip().decode())
""",
    "fixed_param_names": """
out = []
n = c_int(0); E.EFTCAMB_get_num_params(byref(eft), byref(n))
for i in range(n.value):
    buf = create_string_buffer(200)
    E.EFTCAMB_get_param_names(byref(eft), byref(c_int(i + 1)), buf, c_size_t(200))
    out.append(buf.raw.rstrip().decode())
emit(value=out)
""",
    "fixed_param_labels": """
out = []
n = c_int(0); E.EFTCAMB_get_num_params(byref(eft), byref(n))
for i in range(n.value):
    buf = create_string_buffer(200)
    E.EFTCAMB_get_param_labels(byref(eft), byref(c_int(i + 1)), buf, c_size_t(200))
    out.append(buf.raw.rstrip().decode())
emit(value=out)
""",
    "fixed_read_parameters": """
k, v = create_string_buffer(200), create_string_buffer(200)
n = c_int(0)
E.EFTCAMB_get_read_parameters(byref(eft), k, v, byref(n), byref(c_int(0)), c_size_t(200), c_size_t(200))
items = []
for i in range(n.value):
    k, v = create_string_buffer(200), create_string_buffer(200)
    E.EFTCAMB_get_read_parameters(byref(eft), k, v, byref(c_int(i + 1)), byref(c_int(1)), c_size_t(200),
                                  c_size_t(200))
    items.append((k.raw.rstrip().decode(), v.raw.rstrip().decode()))
emit(n=n.value, first=items[:3])
""",
    "fixed_set_model_name": """
buf = create_string_buffer('renamed'.ljust(200).encode(), 200)
E.EFTCAMB_set_model_name(byref(eft), buf, c_size_t(200))
out = create_string_buffer(200)
E.EFTCAMB_get_model_name(byref(eft), out, c_size_t(200))
emit(value=out.raw.rstrip().decode())
""",
}
A0_ONLY = {
    "immutable_bytes_written": """
b = ('{:<200}'.format('')).encode('utf-8')
before = bytes(bytearray(b))  # a real copy (bytes(b) would return b itself)
E.EFTCAMB_get_model_name(byref(eft), b)
emit(bytes_object_mutated=(b != before), before=before[:20].decode(), now=b[:30].decode())
""",
    "multibyte_name_decode": """
eft.set_model_name('Omega' + '\\u03a9' * 150)   # 5 + 2*97 = 199 bytes, byte 200 starts a 2-byte character
try:
    emit(value_len=len(eft.model_name()))
except Exception as e:
    emit(error=type(e).__name__ + ': ' + str(e)[:120])
""",
}

if __name__ == "__main__":
    tree = sys.argv[1] if len(sys.argv) > 1 else "AD"
    steps = dict(STEPS)
    if tree == "A0":
        steps.update(A0_ONLY)
    for name, code in steps.items():
        r = run_snippet(SETUP + code, tree_name=tree, threads=1, timeout=300)
        tail = [t.strip() for t in r["tail"] if "runtime error" in t.lower() or "Error" in t][:2]
        print(f"== {tree} {name}: status={r['status']} {json.dumps(r['result'])[:300]} {tail}", flush=True)
