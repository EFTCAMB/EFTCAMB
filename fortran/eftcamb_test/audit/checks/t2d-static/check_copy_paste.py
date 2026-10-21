"""Copy-paste heuristics over EFTCAMB model / parametrization code.

1. eft_cache%EFT<Name><Order> = ... self%<obj>%<method>(...): the derivative order of <method> must match <Order>
   (V->value, P->first_derivative, PP->second, PPP->third, PPPP->fourth) and the Gamma/Omega index of <obj>
   should match <Name> (e.g. EFTGamma3P from PureEFTGamma3%first_derivative).
2. Parametrization derivative procedures (FirstDerivative ...) whose result is assigned from a sibling
   (e.g. FirstDerivative returning self%...%second_derivative or a variable named for a different order).
3. Consecutive assignments with the same LHS in one procedure (overwritten immediately).
Usage: python check_copy_paste.py
"""

from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fparse  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
ORDER = {"v": 0, "p": 1, "pp": 2, "ppp": 3, "pppp": 4}
METH = {"value": 0, "first_derivative": 1, "second_derivative": 2, "third_derivative": 3, "fourth_derivative": 4}


def main():
    fs = sorted(glob.glob(os.path.join(REPO, "fortran/eftcamb/**/*.f90"), recursive=True))
    fs = [f for f in fs if "02_bundled" not in f]
    n = 0
    for f in fs:
        rel = os.path.relpath(f, REPO)
        for p in fparse.parse_file(f):
            prev = None
            for st in p.body:
                asg = fparse.assignment(st.text)
                if not asg:
                    prev = None
                    continue
                base, lhs, rhs = asg
                rl = rhs.lower()
                m = re.match(r"^eft_cache%eft(omega|gamma(\d)|alpha\w|c|lambda)(v|p|pp|ppp|pppp|dot|dotdot)?$", lhs)
                if m:
                    name = m.group(1)
                    idx = m.group(2)
                    order = m.group(3)
                    calls = re.findall(r"self%(\w+)%(value|first_derivative|second_derivative|third_derivative|fourth_derivative)\b", rl)
                    for obj, meth in calls:
                        msgs = []
                        if order in ORDER and METH[meth] != ORDER[order] and "_ode" not in obj:
                            # derivatives wrt a can legitimately appear in combination; flag only single-term RHS
                            if len(calls) == 1:
                                msgs.append(f"order {order} from {meth}")
                        mi = re.search(r"gamma(\d)", obj)
                        if idx and mi and mi.group(1) != idx:
                            msgs.append(f"index Gamma{idx} from {obj}")
                        if name == "omega" and "omega" not in obj and "gamma" in obj:
                            msgs.append(f"Omega from {obj}")
                        if msgs:
                            n += 1
                            print(f"{rel}:{st.line}: {p.name}: {'; '.join(msgs)} | {st.text[:150]}")
                # consecutive same LHS (immediate overwrite) without use of LHS in RHS
                if prev and prev[0] == lhs and lhs not in rl.replace(" ", ""):
                    n += 1
                    print(f"{rel}:{st.line}: {p.name}: immediate overwrite of {lhs} (prev line {prev[1]}) | {st.text[:120]}")
                prev = (lhs, st.line)
            # parametrization derivative procedures
            pn = p.name.lower()
            want = None
            for k, v in (("fourthderivative", 4), ("thirdderivative", 3), ("secondderivative", 2), ("firstderivative", 1)):
                if pn.endswith(k):
                    want = v
                    break
            if want is not None and p.kind == "function":
                for st in p.body:
                    asg = fparse.assignment(st.text)
                    if asg and asg[0] == pn:
                        for obj, meth in re.findall(r"%(\w+)%(value|first_derivative|second_derivative|third_derivative|fourth_derivative)\b", asg[2].lower()):
                            if METH[meth] != want:
                                n += 1
                                print(f"{rel}:{st.line}: {p.name}: derivative-order mismatch ({meth}) | {st.text[:120]}")
    print("TOTAL", n)


if __name__ == "__main__":
    main()
