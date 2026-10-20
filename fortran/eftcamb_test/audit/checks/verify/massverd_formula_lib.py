"""Helper for massverd_*.py: transcribe the Fortran body of EFTCAMBModelComputeMassFactors into Python statements."""

import re

DEFAULT_SRC = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-"
               "571fb7648a45/scratchpad/audit/trees/A0/fortran/eftcamb/06_abstract_EFTCAMB_model.f90")


def transcribe(path=DEFAULT_SRC, numeric=False):
    text = open(path).read()
    i = text.index("subroutine EFTCAMBModelComputeMassFactors")
    j = text.index("end subroutine EFTCAMBModelComputeMassFactors")
    body = text[i:j]
    body = body[body.index("rho_d = "):]
    joined = ""
    for ln in body.splitlines():
        s = ln.split("!")[0].strip()
        if not s:
            continue
        if joined.endswith("&"):
            joined = joined[:-1] + " " + s
        else:
            joined += "\n" + s
    out = []
    for s in (x for x in joined.split("\n") if x.strip()):
        s = s.replace("eft_cache%", "").replace("eft_par_cache%", "")
        if numeric:
            s = s.replace("_dl", "")
        else:
            s = re.sub(r"(\d+)\._dl", r"\1", s)
            s = re.sub(r"(\d+)\.(?!\d)", r"\1", s)
            s = s.replace("sqrt", "sp.sqrt")
        out.append(s)
    return out
