"""Second pass of the per-k sweep: in per-k reachable procedures list allocate/deallocate statements, pointer
assignments and calls that pass shared objects (self%..., *par_cache*) as actual arguments (possible writes through
intent(out/inout) dummies)."""
import json, os, re, sys, glob, collections
sys.path.insert(0, os.path.dirname(__file__))
import fparse
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
reach = set(json.load(open(sys.argv[1]))["reach"])
fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
for f in fs:
    for p in fparse.parse_file(f):
        if p.name.lower() not in reach or p.in_interface:
            continue
        for st in p.body:
            tl = st.text.lower()
            tag = None
            if re.match(r"^\s*(de)?allocate\s*\(", tl):
                tag = "alloc"
            elif "=>" in tl and not tl.startswith("select") and not tl.startswith("associate"):
                tag = "ptr-assign"
            elif re.match(r"^\s*call\b", tl):
                m = re.match(r"^\s*call\s+([\w%()\s]+?)\s*\((.*)\)\s*$", tl)
                if m:
                    args = fparse.split_top(m.group(2))
                    sh = [a.strip() for a in args if re.match(r"^\s*(self|\w*par_cache|params_cache)\b", a.strip())]
                    callee = re.sub(r"\s", "", m.group(1))
                    if sh and not callee.endswith("%initialize"):
                        tag = "call-shared-arg:" + callee + "(" + ",".join(sh) + ")"
            elif re.search(r"\b(stop|error\s+stop|mpistop)\b", tl):
                tag = "stop"
            if tag:
                print(f"{os.path.relpath(f, REPO)}:{st.line}: {p.name}: {tag[:160]} | {st.text[:120]}")
