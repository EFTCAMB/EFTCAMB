"""Independent re-evaluation (numpy) of the EFTCAMB closed-form quasi-static ratios QS_mu and QS_mu_gen.

The Fortran expressions are read from the frozen A0 source (06_abstract_EFTCAMB_model.f90), translated mechanically to
Python, and split at the single top-level '/' into numerator and denominator, so that one can see whether a NaN is
0/0 (both vanish), x/0 (inf) or something else.  Inputs are arrays of timestep-cache fields (from the harness npz:
eft_<field> background arrays or evol_<field>[k, a] arrays).

Used by verify_qsnan_scan.py.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

SRC = Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0/fortran/eftcamb/06_abstract_EFTCAMB_model.f90"
)


def _sub_body(text: str, name: str) -> str:
    m = re.search(rf"subroutine {name}\(.*?\n(.*?)end subroutine {name}\b", text, re.S)
    if not m:
        raise RuntimeError(name)
    return m.group(1)


def _join(body: str) -> list[str]:
    """Join Fortran continuation lines into logical statements."""
    stmts, cur = [], ""
    for ln in body.splitlines():
        s = ln.split("!")[0].rstrip()
        if not s.strip():
            continue
        cont = s.endswith("&")
        s = s[:-1] if cont else s
        cur += " " + s.strip()
        if not cont:
            stmts.append(cur.strip())
            cur = ""
    return stmts


def _py(expr: str) -> str:
    expr = re.sub(r"eft_par_cache%h0_Mpc", "h0", expr)
    expr = re.sub(r"eft_cache%(\w+)", r"C['\1']", expr)
    expr = re.sub(r"_dl\b", "", expr)
    return expr


def _split_top_div(expr: str) -> tuple[str, str]:
    depth, cut = 0, []
    for i, ch in enumerate(expr):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "/" and depth == 0:
            cut.append(i)
    if len(cut) != 1:
        raise RuntimeError(f"expected one top-level '/', found {len(cut)}")
    return expr[: cut[0]], expr[cut[0] + 1 :]


def load_formulas(src: Path = SRC) -> dict:
    text = src.read_text()
    out = {}
    for name, field in (("EFTCAMBModelComputeQSMu", "QS_mu"), ("EFTCAMBModelComputeQSSigma", "QS_sigma")):
        st = [s for s in _join(_sub_body(text, name)) if s.startswith(f"eft_cache%{field} =")]
        assert len(st) == 1
        rhs = st[0].split("=", 1)[1]
        num, den = _split_top_div(_py(rhs))
        out[field] = (num, den)
    # QS_mu_gen: auxiliary assignments then the ratio
    stm = _join(_sub_body(text, "EFTCAMBModelComputeQSMuSigma"))
    aux = []
    for s in stm:
        m = re.match(r"(\w+) = (.*)", s)
        if m and m.group(1) in ("R0dot", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "M2"):
            aux.append(f"{m.group(1)} = {_py(m.group(2))}")
        if s.startswith("eft_cache%QS_mu_gen ="):
            out["QS_mu_gen"] = _split_top_div(_py(s.split("=", 1)[1]))
    out["QS_mu_gen_aux"] = aux
    return out


def evaluate(F: dict, C: dict, a, k, h0: float) -> dict:
    """Return {field: (num, den)} evaluated with numpy broadcasting."""
    ns = {"C": C, "a": a, "k": k, "h0": h0, "np": np}
    with np.errstate(all="ignore"):
        res = {f: (eval(F[f][0], ns), eval(F[f][1], ns)) for f in ("QS_mu", "QS_sigma")}
        for line in F["QS_mu_gen_aux"]:
            exec(line, ns)
        res["QS_mu_gen"] = (eval(F["QS_mu_gen"][0], ns), eval(F["QS_mu_gen"][1], ns))
    return res
