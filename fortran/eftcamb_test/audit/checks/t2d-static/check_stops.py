"""List stop / error stop / MpiStop sites in EFTCAMB code (and EFTCAMB MOD blocks) with enclosing procedure,
flag those inside per-k reachable procedures (perk.json from check_perk_writes.py)."""
import json, os, re, sys, glob
sys.path.insert(0, os.path.dirname(__file__))
import fparse
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
reach = set(json.load(open(sys.argv[1]))["reach"]) if len(sys.argv) > 1 else set()
fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
n = 0
for f in fs:
    for p in fparse.parse_file(f):
        for st in p.body:
            tl = st.text.lower()
            if re.search(r"(^|\bthen\b|\)\s*)\s*(error\s+stop|stop)\b", tl) or re.search(r"\bcall\s+mpistop\b", tl):
                n += 1
                tag = "PERK" if p.name.lower() in reach else "    "
                print(f"{tag} {os.path.relpath(f, REPO)}:{st.line}: {p.name}: {st.text[:110]}")
print("TOTAL", n)
