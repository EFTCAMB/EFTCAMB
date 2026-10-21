"""Static scan: timestep-cache fields read by model background routines before any caller has filled them.

For every EFTCAMB model implementation of compute_background_EFT_functions, compute_adotoa, compute_H_derivs and
compute_dtauda, list the eft_cache fields the routine reads before writing them itself (textual order, internal
procedures of the routine included when they are called by name), and compare with the fields the callers have set
at that point:

  derivs (equations.f90:2786-2903), full model : bg_EFT -> adotoa -> [nu block] -> H_derivs x2 -> rhoQPQ -> 2nd order
  derivs, designer model                       : adotoa -> [nu block] -> H_derivs x2 -> bg_EFT -> rhoQPQ -> 2nd order
  dtauda (equations.f90:18-50)                 : fresh cache with only a, grhoa2 -> bg_EFT (full) -> compute_dtauda

Usage: python check_cache_read_before_write.py [--json OUT]
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
BINDS = ["compute_background_eft_functions", "compute_adotoa", "compute_h_derivs", "compute_dtauda"]
DERIVS_SET = {"a", "tau", "k", "grhom_t", "grhoa2", "gpresm_t", "grhob_t", "grhoc_t", "grhor_t", "grhog_t"}
NU_SET = {"grhonu_tot", "gpinu_tot", "grhonudot_tot", "gpinudot_tot", "gpresdotm_t"}
NU2_SET = {"gpinudotdot_tot", "gpresdotdotm_t"}
H_SET = {"hdot", "hdotdot", "hdotdotdot"}
BIND_RE = re.compile(r"^\s*procedure\s*(?:,[^:]*)?::\s*(\w+)\s*=>\s*(\w+)", re.I)


def cache_dummies(proc):
    out = set()
    for st in proc.decls:
        if re.search(r"type\s*\(\s*TEFTCAMB_timestep_cache\s*\)", st.text, re.I):
            _, ents = fparse.decl_entities(st.text)
            out |= {e[0].lower() for e in ents}
    return out


def scan_reads(proc, internals, seen=None):
    """Return ordered list of (field, line, kind) with kind 'R' or 'W' for the cache dummies of proc."""
    seen = seen or set()
    names = cache_dummies(proc) or {"eft_cache"}
    events = []
    pat = re.compile(r"\b(" + "|".join(sorted(names)) + r")\s*%\s*(\w+)", re.I)
    for st in proc.body:
        t = st.text
        asg = fparse.assignment(t)
        lhs_field = None
        rhs = t
        if asg:
            _, lhs, rhs = asg
            m = pat.match(lhs.strip())
            if m:
                lhs_field = m.group(2).lower()
            # an 'if (cond) x = y' has the condition in the statement before the lhs
            pre = t[: t.lower().find(lhs.lower())] if lhs.lower() in t.lower() else ""
            rhs = pre + " " + rhs
        for m in pat.finditer(rhs):
            events.append((m.group(2).lower(), st.line, "R"))
        # internal procedures called by name
        for m in re.finditer(r"\b(\w+)\s*\(", t):
            n = m.group(1).lower()
            if n in internals and n not in seen:
                seen.add(n)
                events.extend(scan_reads(internals[n], internals, seen))
        if lhs_field:
            events.append((lhs_field, st.line, "W"))
    return events


def first_reads(events):
    written, out = set(), collections.OrderedDict()
    for f, line, k in events:
        if k == "W":
            written.add(f)
        elif f not in written and f not in out:
            out[f] = line
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json")
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/0[678]*/*.f90")))
    files += sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/06*.f90")))
    res = []
    for f in files:
        binds = {}
        for st in fparse.logical_statements(f):
            m = BIND_RE.match(st.text)
            if m and m.group(1).lower() in BINDS:
                binds[m.group(2).lower()] = m.group(1).lower()
        if not binds:
            continue
        procs = [p for p in fparse.parse_file(f) if not p.in_interface]
        byname = {p.name.lower(): p for p in procs}
        src = open(f).read().lower()
        designer = "extends (eftcamb_designer_model)" in src or "extends(eftcamb_designer_model)" in src or "06p2" in f
        own_writes = {}
        for impl, bname in binds.items():
            p = byname.get(impl)
            if p is not None:
                internals = {q.name.lower(): q for q in procs if (q.parent or "").lower() == impl}
                own_writes[bname] = {x[0] for x in scan_reads(p, internals) if x[2] == "W"}
        for impl, bname in binds.items():
            p = byname.get(impl)
            if p is None:
                continue
            internals = {q.name.lower(): q for q in procs if (q.parent or "").lower() == impl}
            reads = first_reads(scan_reads(p, internals))
            ad = own_writes.get("compute_adotoa", set())
            bg = own_writes.get("compute_background_eft_functions", set())
            hd = own_writes.get("compute_h_derivs", set())
            if bname == "compute_background_eft_functions":
                avail = DERIVS_SET | ((NU_SET | NU2_SET | H_SET | {"adotoa"} | ad | hd) if designer else set())
            elif bname == "compute_adotoa":
                avail = DERIVS_SET | (set() if designer else bg)
            elif bname == "compute_h_derivs":
                avail = DERIVS_SET | NU_SET | {"adotoa"} | ad | (set() if designer else bg)
            else:
                avail = {"a", "grhoa2"} | (set() if designer else bg)
            # fields produced by the model's own EFT-function routine count as available for adotoa/H/dtauda
            suspicious = {k: v for k, v in reads.items() if k not in avail and not k.startswith("eft")}
            res.append(dict(file=os.path.relpath(f, REPO), binding=bname, impl=p.name, line=p.line,
                            reads=list(reads), suspicious=suspicious, designer=designer))
    for r in res:
        if r["suspicious"]:
            s = ", ".join(f"{k}@{v}" for k, v in r["suspicious"].items())
            kind = "designer" if r["designer"] else "full"
            print(f"{r['file']}:{r['line']} [{kind}] {r['binding']} ({r['impl']}): reads before set: {s}")
    print(f"TOTAL routines {len(res)}, with suspicious reads {sum(1 for r in res if r['suspicious'])}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
