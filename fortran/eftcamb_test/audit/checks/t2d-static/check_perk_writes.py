"""Thread-safety sweep: find writes to shared state in EFTCAMB procedures reachable from per-k evolution.

Entry points: the type-bound procedures CAMB calls inside OpenMP parallel loops (equations.f90 derivs/output/
initial/derivst, the IC submodule, CAMBdata_*Arr OMP loops in results.f90).
Reachability: conservative call graph over type-bound bindings (any binding with a given name), direct calls and
function references to known procedures.  Shared state = dummies named self / *par_cache* / params_cache, module
variables, implicit-SAVE locals.

Usage: python check_perk_writes.py [--json OUT] [--verbose]
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))

ENTRY_BINDINGS = [
    "compute_background_eft_functions",
    "compute_secondorder_eft_functions",
    "compute_dtauda",
    "compute_adotoa",
    "compute_h_derivs",
    "compute_rhoqpq",
    "compute_einstein_factors",
    "compute_pi_factors",
    "compute_tensor_factors",
    "compute_stability_factors",
    "compute_mass_factors",
    "compute_qs_mu",
    "compute_qs_sigma",
    "compute_qs_mu_sigma",
    "compute_deltatau",
    "compute_deltacosmict",
    "compute_soundhorizon",
]
ENTRY_PROCS = ["eftcambinitialconditions"]

BIND_RE = re.compile(
    r"^\s*procedure\s*(?:\(\s*\w*\s*\))?\s*(?:,[^:]*)?::\s*(\w+)\s*(?:=>\s*(\w+))?", re.I
)
GENERIC_RE = re.compile(r"^\s*generic\s*(?:,[^:]*)?::\s*(\w+)\s*=>\s*(.*)$", re.I)


def files():
    return sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))


def bindings(fs):
    b = collections.defaultdict(set)
    for f in fs:
        for st in fparse.logical_statements(f):
            m = BIND_RE.match(st.text)
            if m:
                name = m.group(1).lower()
                impl = (m.group(2) or m.group(1)).lower()
                b[name].add(impl)
            m = GENERIC_RE.match(st.text)
            if m:
                for x in m.group(2).split(","):
                    b[m.group(1).lower()].add(("@" + x.strip().lower()))
    # resolve generic -> specific bindings
    out = collections.defaultdict(set)
    for k, v in b.items():
        for x in v:
            if x.startswith("@"):
                out[k] |= {y for y in b.get(x[1:], set()) if not y.startswith("@")}
            else:
                out[k].add(x)
    return out


def calls_of(proc, binds, procnames):
    """Return set of callee procedure names (lowercase) referenced in body."""
    out = set()
    for st in proc.body:
        t = st.text
        tl = t.lower()
        # method references a%b%name(
        for m in re.finditer(r"%\s*(\w+)\s*\(", tl):
            n = m.group(1)
            if n in binds:
                out |= binds[n]
        for m in re.finditer(r"%\s*(\w+)\s*$", tl):  # call x%name without parens
            n = m.group(1)
            if n in binds and tl.lstrip().startswith("call"):
                out |= binds[n]
        m = re.match(r"^\s*call\s+(\w+)\s*(\(|$)", tl)
        if m:
            out.add(m.group(1))
        for m in re.finditer(r"(?<![%\w])(\w+)\s*\(", tl):
            n = m.group(1)
            if n in procnames:
                out.add(n)
        # procedure arguments passed to integrators / root finders (function names as actual args)
        for m in re.finditer(r"[(,]\s*(\w+)\s*(?=[,)])", tl):
            n = m.group(1)
            if n in procnames:
                out.add(n)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    fs = files()
    procs = []
    for f in fs:
        procs.extend(p for p in fparse.parse_file(f) if not p.in_interface)
    byname = collections.defaultdict(list)
    for p in procs:
        byname[p.name.lower()].append(p)
    binds = bindings(fs)
    procnames = set(byname)
    graph = {}
    for name, ps in byname.items():
        s = set()
        for p in ps:
            s |= calls_of(p, binds, procnames)
        graph[name] = {x for x in s if x in procnames and x != name}
    # reachability
    start = set()
    for b in ENTRY_BINDINGS:
        start |= binds.get(b, set())
    start |= set(ENTRY_PROCS)
    start &= procnames
    reach = set()
    parent = {}
    stack = list(start)
    for s in start:
        parent[s] = None
    while stack:
        n = stack.pop()
        if n in reach:
            continue
        reach.add(n)
        for c in graph.get(n, ()):
            if c not in reach:
                parent.setdefault(c, n)
                stack.append(c)
    print(f"procedures: {len(procs)}, per-k reachable: {len(reach)} (entries {len(start)})")
    # module variables
    modvars = {}
    sys.path.insert(0, os.path.dirname(__file__))
    from check_module_vars import module_level_decls  # noqa: E402

    for f in fs:
        for line, n, spec, v in module_level_decls(f):
            modvars[n] = f"{os.path.relpath(f, REPO)}:{line}"
    rows = []
    for name in sorted(reach):
        for p in byname[name]:
            dummies = set(p.args)
            locals_ = set()
            saved = set()
            for d in p.decls:
                if fparse.TYPE_RE.match(d.text.lower()):
                    spec, ents = fparse.decl_entities(d.text)
                    for n, v in ents:
                        if n not in dummies:
                            locals_.add(n)
                            if v is not None and "parameter" not in spec.lower():
                                saved.add(n)
            for st in p.body:
                asg = fparse.assignment(st.text)
                if not asg:
                    continue
                base, lhs, rhs = asg
                kind = None
                if base == "self" or "par_cache" in base or base == "params_cache":
                    if base in dummies:
                        kind = "shared-dummy"
                elif base in saved:
                    kind = "implicit-save"
                elif base in modvars and base not in locals_ and base not in dummies:
                    kind = "module-var"
                elif base not in locals_ and base not in dummies and p.parent is None and base != (p.result or ""):
                    kind = "undeclared?"
                if kind:
                    chain = []
                    n = name
                    while n is not None and len(chain) < 8:
                        chain.append(n)
                        n = parent.get(n)
                    rows.append(
                        dict(
                            file=os.path.relpath(p.file, REPO),
                            line=st.line,
                            proc=p.name,
                            kind=kind,
                            lhs=lhs,
                            chain=" <- ".join(chain),
                        )
                    )
    for r in rows:
        if r["kind"] == "undeclared?" and not a.verbose:
            continue
        print(f"{r['file']}:{r['line']}: [{r['kind']}] {r['proc']}: {r['lhs']} = ...   via {r['chain']}")
    if a.verbose:
        for n in sorted(reach):
            print("REACH", n)
    if a.json:
        json.dump(dict(rows=rows, reach=sorted(reach)), open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
