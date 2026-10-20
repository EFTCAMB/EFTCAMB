"""t8b check: equations.f90 dtauda (EFTCAMB MOD) evaluates sqrt(3/grhoa2) before testing grhoa2 <= 0.
In a build with -ffpe-trap=invalid (AD) a GR model whose expansion rate becomes imaginary (loitering closed model)
then dies with SIGFPE instead of raising CAMBError('Universe stops expanding ...'). Usage: python dtauda_recollapse.py"""

import sys


def main(args):
    import camb

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.0, ombh2=0.022, omch2=0.12, omk=-1.8, mnu=0.0, num_massive_neutrinos=0)
    try:
        camb.get_background(pars, no_thermo=True)
        return {"result": "no error"}
    except Exception as e:  # noqa: BLE001
        return {"result": repr(e)[:160]}


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    for t in sys.argv[1:] or ["A0", "AD", "B0"]:
        r = run_in_tree(t, HERE / "dtauda_recollapse.py", {}, threads=1)
        print(t, json.dumps(r)[:600])
