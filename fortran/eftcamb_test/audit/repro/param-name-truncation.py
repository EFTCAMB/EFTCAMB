"""Repro (t4a-param): ParametrizedFunction1DInitFromFile (04_abstract_parametrizations_1D.f90:279-296) reads every
parameter with a key stored in character(len=EFT_names_max_length) = 20, so parameter names longer than 20
characters are silently truncated and the documented ini key is never read (default 0 is used instead).
Affected in the code: the 3+1 designer functions 'Designer3p1sigma1_ODE' (21 chars), 'Designer3p1rho1N_ODE' (20),
'Designer3p1rho2_ODE' (19, with 'Exp'/'_v1'... suffixes), and RPH '*_ODE' functions with steplog/Fourier suffixes.

With --model also runs RPH 3_RPH_KT_OmegaDE_1.ini with the kineticity given as a steplog function of Omega_DE
(RPHkineticitymodel_ODE = 8, v1 = 0.01, v2 = 0.2, at = 0.3, RPHkineticity_ODE_delta = 0.5): the 23-character
documented key RPHkineticity_ODE_delta is never read (delta = 0 -> NaN); then the same run with the truncated
key RPHkineticity_ODE_de = 0.5, which is what the code actually reads.  Expected: both ok and identical.

usage: /Users/marco/Software/miniforge3/bin/python param-name-truncation.py [--model]
"""

import subprocess
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "t4a-param" / "repro" / Path(__file__).stem
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program repro
use precision
use IniObjects
use EFT_def
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i, flag
real(dl) :: v
character(len=200) :: longname
character(len=EFT_names_max_length) :: key
character(len=512) :: path, fname, sflag
call get_command_argument(1, path)
call get_command_argument(2, fname)
call get_command_argument(3, sflag)
read(sflag, *) flag
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, flag, trim(fname), 'x', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
do i = 1, f%parameter_number
    call f%parameter_names(i, longname)
    call f%parameter_names(i, key)
    call f%parameter_value(i, v)
    write(*, '(a,"|",a,"|",es14.6)') trim(longname), trim(key), v
end do
end program
"""


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    for fname, flag, keys in [
        ("EFTOmega", 3, {"EFTOmega0": 0.7, "EFTOmegaExp": 2.0}),  # control
        ("Designer3p1sigma1_ODE", 3, {"Designer3p1sigma1_ODE0": 0.7, "Designer3p1sigma1_ODEExp": 2.0}),
        ("Designer3p1rho2_ODE", 3, {"Designer3p1rho2_ODE0": 0.7, "Designer3p1rho2_ODEExp": 2.0}),
        ("RPHbraiding_ODE", 8, {"RPHbraiding_ODE_v1": 0.1, "RPHbraiding_ODE_v2": 0.2, "RPHbraiding_ODE_at": 0.3,
                                "RPHbraiding_ODE_delta": 0.4}),
        ("RPHkineticity_ODE", 8, {"RPHkineticity_ODE_v1": 0.1, "RPHkineticity_ODE_v2": 0.2,
                                  "RPHkineticity_ODE_at": 0.3, "RPHkineticity_ODE_delta": 0.4}),
    ]:
        (WORK / "r.ini").write_text("".join(f"{k} = {v}\n" for k, v in keys.items()))
        out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini"), fname, str(flag)], capture_output=True,
                             text=True, check=True).stdout
        print(f"function name {fname!r} ({len(fname)} chars), allocator flag {flag}:")
        for line in out.splitlines():
            name, key, v = line.split("|")
            exp = keys.get(name)
            mark = "" if exp is not None and abs(float(v) - exp) < 1e-12 else "   <-- ini value not read"
            print(f"   param {name:28s} key used {key!r:24s} read {float(v):8.3f}  expected {exp}{mark}")


def model():
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
    import harness
    ini = harness.PARAMS_DIR / "3_RPH_KT_OmegaDE_1.ini"
    base = {"RPHkineticitymodel_ODE": "8", "RPHkineticity_ODE_v1": "0.01", "RPHkineticity_ODE_v2": "0.2",
            "RPHkineticity_ODE_at": "0.3"}
    runs = {"documented key RPHkineticity_ODE_delta = 0.5": {**base, "RPHkineticity_ODE_delta": "0.5"},
            "truncated key RPHkineticity_ODE_de = 0.5": {**base, "RPHkineticity_ODE_de": "0.5"}}
    files = {}
    for i, (name, ov) in enumerate(runs.items()):
        f = WORK / f"model{i}.npz"
        st = harness.run_model("A0", ini, f, overrides=ov, threads=1)
        tail = [t.strip() for t in (st.get("log_tail") or []) if "NaN" in t or "unstable" in t][:2]
        files[name] = f if st["status"] == "ok" else None
        print(f"{name:46s} status {st['status']}  {' / '.join(tail) if st['status'] != 'ok' else ''}")
    a, b = list(files.values())
    if a and b:
        c = harness.compare(harness.load(a), harness.load(b))
        print(f"difference between the two runs: worst Cl rel diff = {c['worst_cl']:.3e} (expected 0)")


if __name__ == "__main__":
    import sys
    main()
    if "--model" in sys.argv:
        model()
