"""Minimal free-form Fortran parser for the t2d-static audit sweep.

Splits files into logical statements (comments stripped, continuations joined, ';' split),
groups them into program units (module / submodule / subroutine / function, nested contains)
and records, for each procedure: dummy arguments, local declarations (with initialisers),
executable statements, assignments (LHS base name) and calls.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

TYPE_RE = re.compile(
    r"^\s*(real|integer|logical|character|complex|double\s+precision|type\s*\(|class\s*\(|procedure\s*\()", re.I
)
PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure|module|non_recursive)\s+)*"
    r"(?:(?:real|integer|logical|character|complex|double\s+precision)(?:\s*\([^)]*\))?\s+|type\s*\([^)]*\)\s+)?"
    r"(?:(?:pure|elemental|recursive|impure|module)\s+)*"
    r"(subroutine|function)\s+(\w+)\s*(\(([^)]*)\))?(.*)$",
    re.I,
)
END_RE = re.compile(r"^\s*end\s*(subroutine|function|module|submodule|program|procedure)\b(\s+\w+)?\s*$", re.I)
END_PLAIN_RE = re.compile(r"^\s*end\s*$", re.I)
MODULE_RE = re.compile(r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)(\w+)\s*$", re.I)
SUBMODULE_RE = re.compile(r"^\s*submodule\s*\(([^)]*)\)\s*(\w+)", re.I)
TYPE_DEF_RE = re.compile(r"^\s*type\s*(?:,[^:]*)?(?:::)?\s*(\w+)\s*$", re.I)
INTERFACE_RE = re.compile(r"^\s*(abstract\s+)?interface\b", re.I)
MODPROC_RE = re.compile(r"^\s*module\s+procedure\s+(\w+)\s*$", re.I)


@dataclass
class Stmt:
    file: str
    line: int  # first physical line
    text: str  # logical statement, comments removed
    raw: str  # raw text of the physical lines joined


@dataclass
class Proc:
    name: str
    kind: str
    file: str
    line: int
    module: str
    args: list[str]
    result: str | None
    decls: list[Stmt] = field(default_factory=list)
    body: list[Stmt] = field(default_factory=list)
    parent: str | None = None
    end_line: int = 0
    in_interface: bool = False


def strip_comment(line: str) -> str:
    out = []
    q = None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        else:
            if ch in "'\"":
                q = ch
                out.append(ch)
            elif ch == "!":
                break
            else:
                out.append(ch)
    return "".join(out)


def logical_statements(path: str) -> list[Stmt]:
    lines = Path(path).read_text(errors="replace").split("\n")
    stmts: list[Stmt] = []
    buf = ""
    raw = ""
    start = 0
    for i, ln in enumerate(lines, 1):
        s = ln.strip()
        if s.startswith("#"):
            # preprocessor line: keep as its own statement
            if not buf:
                stmts.append(Stmt(path, i, s, ln))
            continue
        c = strip_comment(ln).rstrip()
        if not buf:
            start = i
            raw = ""
        raw += ln + "\n"
        cs = c.strip()
        if cs.startswith("&"):
            cs = cs[1:]
        if cs.endswith("&"):
            buf += cs[:-1] + " "
            continue
        buf += cs
        if buf.strip():
            # split on ; outside strings
            parts = []
            q = None
            cur = ""
            for ch in buf:
                if q:
                    cur += ch
                    if ch == q:
                        q = None
                elif ch in "'\"":
                    q = ch
                    cur += ch
                elif ch == ";":
                    parts.append(cur)
                    cur = ""
                else:
                    cur += ch
            parts.append(cur)
            for p in parts:
                if p.strip():
                    stmts.append(Stmt(path, start, p.strip(), raw))
        buf = ""
    return stmts


def is_decl(text: str) -> bool:
    t = text.strip().lower()
    # an assignment to a variable whose name happens to be an attribute keyword (value = ..., data = ...)
    if re.match(r"^[a-z_]\w*\s*(\([^=]*\))?\s*(%\s*\w+\s*(\([^=]*\))?\s*)*=(?![=>])", t) and "::" not in t:
        return False
    if TYPE_RE.match(t):
        # avoid "real(dl) function" etc
        if re.search(r"\bfunction\b", t) and "::" not in t:
            return False
        return True
    return bool(
        re.match(
            r"^(implicit|use|save|data|parameter|external|intrinsic|dimension|allocatable|pointer|target|optional|"
            r"import|equivalence|common|namelist|private|public|intent|procedure\s*::|interface|abstract interface|"
            r"enum|enumerator|end enum|contiguous|asynchronous|volatile|protected|value|format)\b",
            t,
        )
    )


def parse_file(path: str) -> list[Proc]:
    stmts = logical_statements(path)
    procs: list[Proc] = []
    stack: list[Proc] = []
    module = ""
    in_type = 0
    iface_depth = 0
    for st in stmts:
        t = st.text
        tl = t.lower().strip()
        if tl.startswith("#"):
            if stack:
                stack[-1].body.append(st)
            continue
        m = MODULE_RE.match(t)
        if m and not stack:
            module = m.group(1)
            continue
        m = SUBMODULE_RE.match(t)
        if m and not stack:
            module = m.group(2)
            continue
        if INTERFACE_RE.match(tl) and not tl.startswith("interface operator") and "assignment" not in tl:
            iface_depth += 1
            continue
        if re.match(r"^\s*end\s+interface\b", tl):
            iface_depth -= 1
            continue
        if in_type:
            if re.match(r"^\s*end\s*type\b", tl):
                in_type -= 1
            continue
        if TYPE_DEF_RE.match(t) and not tl.startswith("type(") and not re.match(r"^type\s+is\b", tl):
            # derived type definition
            if "::" in t or re.match(r"^\s*type\s+\w+\s*$", t, re.I) or re.match(r"^\s*type\s*,", t, re.I):
                in_type += 1
                continue
        m = PROC_START_RE.match(t)
        if m and not tl.startswith("end") and "::" not in t.split("(")[0]:
            kind = m.group(1).lower()
            name = m.group(2)
            args = [a.strip().lower() for a in (m.group(4) or "").split(",") if a.strip()]
            rest = m.group(5) or ""
            rm = re.search(r"result\s*\(\s*(\w+)\s*\)", rest, re.I)
            result = rm.group(1).lower() if rm else (name.lower() if kind == "function" else None)
            p = Proc(name, kind, path, st.line, module, args, result, parent=stack[-1].name if stack else None)
            p.in_interface = iface_depth > 0
            procs.append(p)
            stack.append(p)
            p._in_decl = True
            continue
        m = MODPROC_RE.match(t)
        if m and not stack and iface_depth == 0:
            # separate module procedure (submodule implementation)
            p = Proc(m.group(1), "subroutine", path, st.line, module, [], None)
            procs.append(p)
            stack.append(p)
            p._in_decl = True
            continue
        if END_RE.match(t) or END_PLAIN_RE.match(t):
            kw = (END_RE.match(t).group(1) or "").lower() if END_RE.match(t) else ""
            if kw in ("module", "submodule", "program"):
                module = ""
                continue
            if stack and kw in ("subroutine", "function", "procedure", ""):
                p = stack.pop()
                p.end_line = st.line
            continue
        if not stack:
            continue
        p = stack[-1]
        if tl == "contains":
            p._in_decl = False
            p._contains = True
            continue
        if getattr(p, "_contains", False):
            continue
        if p._in_decl and is_decl(t):
            p.decls.append(st)
        else:
            p._in_decl = False
            p.body.append(st)
    return procs


def split_top(s: str, sep: str = ",") -> list[str]:
    out, cur, depth, q = [], "", 0, None
    for ch in s:
        if q:
            cur += ch
            if ch == q:
                q = None
            continue
        if ch in "'\"":
            q = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == sep and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    return out


def decl_entities(text: str) -> tuple[str, list[tuple[str, str | None]]]:
    """Return (type-spec+attrs, [(name, initialiser or None)]) for a type declaration statement."""
    if "::" in text:
        spec, ents = text.split("::", 1)
    else:
        m = re.match(r"^\s*((?:real|integer|logical|character|complex|double\s+precision)(?:\s*\([^)]*\))?(?:\*\d+)?)\s+(.*)$", text, re.I)
        if not m:
            return text, []
        spec, ents = m.group(1), m.group(2)
    res = []
    for e in split_top(ents):
        e = e.strip()
        if not e:
            continue
        if "=>" in e:
            n, v = e.split("=>", 1)
            res.append((re.match(r"\s*(\w+)", n).group(1).lower(), "=>" + v.strip()))
        elif "=" in e:
            # careful: '=' inside parentheses (kind=) only in spec, entity is name[(dims)] [= init]
            parts = split_top(e, "=")
            n = parts[0]
            v = "=".join(parts[1:])
            mm = re.match(r"\s*(\w+)", n)
            if mm:
                res.append((mm.group(1).lower(), v.strip()))
        else:
            mm = re.match(r"\s*(\w+)", e)
            if mm:
                res.append((mm.group(1).lower(), None))
    return spec, res


ASSIGN_RE = re.compile(r"^\s*(?:\w+\s*:\s*)?(?:if\s*\((?P<cond>.*)\)\s*)?(?P<lhs>[A-Za-z_]\w*(?:\s*\([^=]*?\))?(?:\s*%\s*\w+(?:\s*\([^=]*?\))?)*)\s*=(?![=>])(?P<rhs>.*)$")


def assignment(text: str) -> tuple[str, str, str] | None:
    """Return (base_name, full_lhs, rhs) for an assignment statement, else None."""
    t = text.strip()
    tl = t.lower()
    if re.match(r"^(do|if\s*\(.*\)\s*then|else|select|case|where|forall|call|print|write|read|open|close|allocate|deallocate|return|stop|error\s+stop|exit|cycle)\b", tl):
        # single-line 'if (cond) x = y' handled below
        m = re.match(r"^if\s*\(", tl)
        if not m or tl.endswith("then"):
            if not tl.startswith("do ") or "=" not in tl:
                return None
            return None
        # find matching paren for if condition
        depth = 0
        for i, ch in enumerate(t):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    rest = t[i + 1:].strip()
                    return assignment(rest) if rest and not rest.lower().startswith("then") else None
        return None
    if is_decl(t):
        return None
    # LHS: up to first '=' at depth 0 that is not '==', '=>', '<=', '>=', '/='
    depth = 0
    q = None
    for i, ch in enumerate(t):
        if q:
            if ch == q:
                q = None
            continue
        if ch in "'\"":
            q = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif ch == "=" and depth == 0:
            prev = t[i - 1] if i else ""
            nxt = t[i + 1] if i + 1 < len(t) else ""
            if nxt in "=>" or prev in "<>=/":
                return None
            lhs = t[:i].strip()
            rhs = t[i + 1:].strip()
            m = re.match(r"^([A-Za-z_]\w*)", lhs)
            if not m:
                return None
            return m.group(1).lower(), re.sub(r"\s+", "", lhs.lower()), rhs
    return None


def all_procs(paths: list[str]) -> list[Proc]:
    out = []
    for p in paths:
        out.extend(parse_file(p))
    return out
