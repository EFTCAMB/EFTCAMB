"""Check the bind(c) Fortran interfaces in 008p0_Horndeski.f90 against the C prototypes in hdsk_coefficients*.c.

For each C function: argument order, count, names, scalar-by-value (double) vs by-reference (double *), and that every
Fortran dummy is real(c_double) with matching value attribute. Also checks each Fortran call site passes the right
number of actual arguments.

Usage: python check_hdsk_interfaces.py [REPO]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[5]
D = REPO / "fortran" / "eftcamb" / "08f_full_models"


def c_protos() -> dict[str, list[tuple[str, str]]]:
    res = {}
    for f in ("hdsk_coefficients.c", "hdsk_coefficients_f.c"):
        src = (D / f).read_text()
        for m in re.finditer(r"^(\w+)\s+(\w+)\s*\(([^)]*)\)\s*\{", src, re.M):
            ret, name, args = m.groups()
            lst = []
            for a in args.split(","):
                a = " ".join(a.split())
                mm = re.match(r"(const\s+)?(\w+)\s*(\*?)\s*(\w+)$", a)
                assert mm, (name, a)
                lst.append((mm.group(4), mm.group(2) + mm.group(3)))
            res[name] = (ret, lst, f)
    return res


def f_interfaces(src: str) -> dict[str, list[tuple[str, str, bool]]]:
    block = src[src.index("    interface") : src.index("end interface")]
    res = {}
    for m in re.finditer(r"subroutine\s+(\w+)\s*\((.*?)\)\s*bind\s*\(\s*c\s*(?:,\s*name\s*=\s*['\"](\w+)['\"])?\s*\)(.*?)end subroutine",
                         block, re.S | re.I):
        name, args, cname, body = m.groups()
        args = [a.strip() for a in args.replace("&", " ").replace("\n", " ").split(",") if a.strip()]
        decl = {}
        for line in body.splitlines():
            line = line.split("!", 1)[0]
            mm = re.match(r"\s*(real|integer|logical|type)\s*\((\w+)\)\s*(,\s*[^:]*)?::\s*(.*)$", line, re.I)
            if mm:
                typ, kind, attrs, names = mm.groups()
                for n in names.split(","):
                    n = n.strip()
                    dims = ""
                    if "(" in n:
                        n, dims = n.split("(", 1)
                    decl[n.strip().lower()] = (f"{typ}({kind})", (attrs or "").replace(" ", "").lower(), dims)
        res[(cname or name).lower()] = [(a.lower(), *decl.get(a.lower(), ("UNDECLARED", "", ""))) for a in args]
    return res


def call_sites(src: str, names: list[str]) -> dict[str, list[int]]:
    """Return number of actual args at each call of each interface routine."""
    out = {}
    flat = re.sub(r"&\s*\n\s*&?", " ", src)
    for n in names:
        counts = []
        for m in re.finditer(rf"call\s+{n}\s*\(", flat, re.I):
            i, depth = m.end(), 1
            j = i
            while depth:
                ch = flat[j]
                depth += ch == "("
                depth -= ch == ")"
                j += 1
            inner = flat[i : j - 1]
            depth, parts, cur = 0, [], ""
            for ch in inner:
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                if ch == "," and depth == 0:
                    parts.append(cur)
                    cur = ""
                else:
                    cur += ch
            parts.append(cur)
            counts.append(len(parts))
        out[n] = counts
    return out


def main() -> int:
    cp = c_protos()
    fsrc = (D / "008p0_Horndeski.f90").read_text()
    fi = f_interfaces(fsrc)
    problems = 0
    for name, (ret, cargs, cfile) in cp.items():
        if name not in fi:
            print(f"[MISSING] C function {name} ({cfile}) has no Fortran interface")
            problems += 1
            continue
        fargs = fi[name]
        if ret != "void":
            print(f"[RET] {name} returns {ret} but Fortran declares subroutine")
            problems += 1
        if len(cargs) != len(fargs):
            print(f"[COUNT] {name}: C {len(cargs)} args, Fortran {len(fargs)}")
            problems += 1
        for i, ((cn, ct), (fn, ft, fattr, fdims)) in enumerate(zip(cargs, fargs)):
            issues = []
            if cn.lower() != fn:
                issues.append(f"name C={cn} F={fn}")
            if ct not in ("double", "double*"):
                issues.append(f"C type {ct}")
            if ft != "real(c_double)":
                issues.append(f"F type {ft}")
            byval = "value" in fattr
            if (ct == "double") != byval:
                issues.append(f"value mismatch C={ct} F_value={byval}")
            if fdims:
                issues.append(f"F array dims ({fdims}")
            if "intent(in)" in fattr and ct.endswith("*"):
                issues.append("output declared intent(in)")
            if issues:
                problems += 1
                print(f"[ARG] {name} #{i + 1}: " + "; ".join(issues))
        print(f"{name:42s} C args={len(cargs):3d} (out={sum(t.endswith('*') for _, t in cargs)}) F args={len(fargs):3d}")
    for name in fi:
        if name not in cp:
            print(f"[EXTRA] Fortran interface {name} has no C definition")
            problems += 1
    # call sites across EFTCAMB
    allsrc = "\n".join(p.read_text() for p in (REPO / "fortran" / "eftcamb").rglob("*.f90"))
    cs = call_sites(allsrc, list(fi))
    for n, counts in cs.items():
        bad = [c for c in counts if c != len(fi[n])]
        print(f"call {n:42s} sites={len(counts)} nargs={sorted(set(counts))} expected={len(fi[n])}"
              + ("  MISMATCH" if bad else ""))
        problems += bool(bad)
    print(f"problems={problems}")
    return problems


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
