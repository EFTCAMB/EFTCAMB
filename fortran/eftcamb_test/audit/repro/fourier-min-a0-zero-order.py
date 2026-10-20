"""Repro (t4a-param): Fourier parametrization (allocator flag 7) with EFT*_min_a0 = T and both orders 0
(EFT*_Fourier_order_cos = EFT*_Fourier_order_sin = 0, accepted by init_func_from_file).
FourierMinimum (04p011_fourier_parametrizations_1D.f90:319-336) samples num_points = 8*max(order_cos, order_sin)+1 = 1
points, so the root-search loop never runs and the minimum stays at its 1e16 initialisation; the offset becomes
-1e16 and value(a) = a0 - 1e16.  Expected: the function is the constant a0 (its own minimum).

usage: /Users/marco/Software/miniforge3/bin/python fourier-min-a0-zero-order.py
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
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err
character(len=512) :: path
call get_command_argument(1, path)
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 7, 'EFTOmega', '\Omega', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
write(*, '(3es26.17)') f%value(0.3_dl), f%value(0.7_dl), f%value(1._dl)
end program
"""


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    for min_a0 in ("F", "T"):
        ini = (f"EFTOmega_Fourier_order_cos = 0\nEFTOmega_Fourier_order_sin = 0\nEFTOmega_min_a0 = {min_a0}\n"
               "EFTOmega_a0 = 0.2\nEFTOmega_phase = 0.0\n")
        (WORK / "r.ini").write_text(ini)
        out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini")], capture_output=True, text=True, check=True)
        vals = list(map(float, out.stdout.split()))
        print(f"min_a0 = {min_a0}: value(0.3, 0.7, 1.0) observed = {vals}   expected = [0.2, 0.2, 0.2]")


if __name__ == "__main__":
    main()
