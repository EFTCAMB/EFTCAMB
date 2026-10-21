"""Token-normalised diff of a Fortran file between two git refs (ignores whitespace, case, line continuations)."""
import difflib
import re
import subprocess
import sys

ref_a, ref_b, path = sys.argv[1:4]


def norm(ref):
    txt = subprocess.run(["git", "show", f"{ref}:{path}"], capture_output=True, text=True, check=True).stdout
    out = []
    for line in txt.splitlines():
        code = re.sub(r"!.*$", "", line) if not line.lstrip().startswith("!$") else line
        out.append(code)
    s = " ".join(out).replace("&", " ")
    s = re.sub(r"\s+", "", s).lower()
    s = s.replace("real(", "real(")
    # split into statements roughly at keywords for readability
    return re.split(r"(?=\b(?:if|call|end|do|else|write|return)\b)", s)


a, b = norm(ref_a), norm(ref_b)
for l in difflib.unified_diff(a, b, lineterm="", n=0):
    print(l[:300])
