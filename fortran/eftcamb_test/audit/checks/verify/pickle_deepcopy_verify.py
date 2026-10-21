"""Independent verification of pickle-deepcopy-eftcamb-segfault.

Cases (each in its own subprocess, 1 CPU slot):
  gr_pickle      : EFTflag=0 params pickled -> must work (control: generic pickling is fine)
  eft_state      : inspect __getstate__ keys/EFTCAMB entry of an EFT params (what gets serialised)
  eft_restore    : unpickle, then re-run EFTCAMB.initialize_parameters with the original read_parameters()
                   -> if root cause is the missing model, H(z=1) should match the original
  eft_numparams  : unpickle and call EFTCAMB.num_params() (claimed to crash too)
  eft_results    : unpickle and call camb.get_results (another entry point)
  eft_repr       : eval(repr(cp)) round trip
Usage: python pickle_deepcopy_verify.py <tree>
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")


def mk(camb, eft=True):
    kw = dict(H0=70, ombh2=0.0226, omch2=0.112)
    if eft:
        kw.update(feedback_level=0, EFTflag=1, PureEFTmodel=1, PureEFTmodelOmega=1, EFTwDE=0, EFTOmega0=1.0)
    return camb.set_params(**kw)


def child(tree, mode):
    sys.path.insert(0, tree)
    import pickle

    import camb

    if mode == "gr_pickle":
        cp = mk(camb, eft=False)
        h0 = camb.get_background(cp).hubble_parameter(1.0)
        c2 = pickle.loads(pickle.dumps(cp))
        print(f"GR orig H1={h0:.12f} pickled H1={camb.get_background(c2).hubble_parameter(1.0):.12f}")
    elif mode == "eft_state":
        cp = mk(camb)
        st = cp.__getstate__()
        print("EFTCAMB in state:", "EFTCAMB" in st, type(st.get("EFTCAMB")).__name__)
        e = st.get("EFTCAMB")
        if e is not None:
            print(" EFTCAMB state keys:", sorted(e.__getstate__().keys())[:5], "... n=", len(e.__getstate__()))
    elif mode == "eft_restore":
        cp = mk(camb)
        h0 = camb.get_background(cp).hubble_parameter(1.0)
        pars = cp.EFTCAMB.read_parameters()
        c2 = pickle.loads(pickle.dumps(cp))
        c2.EFTCAMB.initialize_parameters(c2, pars, print_header=False)
        h2 = camb.get_background(c2).hubble_parameter(1.0)
        print(f"EFT orig H1={h0:.12f} restored+reinit H1={h2:.12f} diff={h2 - h0:.3e}")
    elif mode == "eft_numparams":
        cp = mk(camb)
        c2 = pickle.loads(pickle.dumps(cp))
        print("num_params orig", cp.EFTCAMB.num_params(), flush=True)
        print("num_params pickled", c2.EFTCAMB.num_params(), flush=True)
    elif mode == "eft_results":
        cp = mk(camb)
        c2 = pickle.loads(pickle.dumps(cp))
        c2.set_for_lmax(500)
        r = camb.get_results(c2)
        print("get_results ok", r.get_cmb_power_spectra()["total"][100, 0])
    elif mode == "eft_repr":
        cp = mk(camb)
        s = repr(cp)
        c2 = eval(s, {"camb": camb})
        print("repr EFTflag", c2.EFTCAMB.EFTflag, flush=True)
        print("repr H1", camb.get_background(c2).hubble_parameter(1.0), flush=True)


def main():
    import harness

    tree = sys.argv[1]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    print("tree", tree)
    for mode in ("gr_pickle", "eft_state", "eft_restore", "eft_numparams", "eft_results", "eft_repr"):
        with harness.slot(1):
            p = subprocess.run([sys.executable, __file__, "_child", tree, mode], env=env, capture_output=True,
                               text=True, timeout=600)
        st = "ok" if p.returncode == 0 else (signal.Signals(-p.returncode).name if p.returncode < 0
                                            else f"rc={p.returncode}")
        print(f"[{mode}] status={st}\n  " + p.stdout.strip().replace("\n", "\n  "))
        if p.returncode != 0:
            print("  stderr tail:", p.stderr.strip()[-600:].replace("\n", "\n    "))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_child":
        child(sys.argv[2], sys.argv[3])
    else:
        main()
