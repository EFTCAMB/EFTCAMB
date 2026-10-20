"""Independent check of eft-getters-gr-segfault: different trigger paths, several trees, backtraces."""
import os, signal, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness

S = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREES = {"A0": f"{S}/trees/A0", "AD": f"{S}/trees/AD", "B1": f"{S}/trees/B1", "A0fix": f"{S}/work/t9-wrapper/A0fix"}
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
case = sys.argv[2]
p = camb.CAMBparams()
p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
p.set_for_lmax(300, lens_potential_accuracy=0)
if case == 'num_params_gr':
    print('VAL', p.EFTCAMB.num_params(), flush=True)
elif case == 'num_params_after_valid_then_gr':
    # valid EFT init (Pure EFT, Lambda-like), then re-init as GR: reset in place should drop the model
    p.EFTCAMB.initialize_parameters(p, {'EFTflag': 1, 'PureEFTmodel': 1, 'EFTwDE': 0, 'PureEFTmodelOmega': 0}, print_header=False)
    print('after valid: n =', p.EFTCAMB.num_params(), flush=True)
    p.EFTCAMB.initialize_parameters(p, {'EFTflag': 0}, print_header=False)
    print('VAL', p.EFTCAMB.num_params(), flush=True)
elif case == 'setflag_get_results':
    # different path: no initialize_parameters at all, EFTflag poked through the ctypes field
    p.EFTCAMB.EFTflag = 1
    try:
        r = camb.get_results(p); print('VAL results', r.tau0, flush=True)
    except Exception as e:
        print('VAL exception', repr(e)[:120], flush=True)
elif case == 'failed_init_get_results':
    try:
        p.EFTCAMB.initialize_parameters(p, {'EFTflag': 1, 'PureEFTmodel': 9}, print_header=False)
    except ValueError as e:
        print('raised', e, flush=True)
    print('EFTflag', p.EFTCAMB.EFTflag, flush=True)
    try:
        r = camb.get_results(p); print('VAL results', r.tau0, flush=True)
    except Exception as e:
        print('VAL exception', repr(e)[:120], flush=True)
"""
cases = ["num_params_gr", "num_params_after_valid_then_gr", "setflag_get_results", "failed_init_get_results"]
env = dict(os.environ, OMP_NUM_THREADS="1"); env.pop("PYTHONPATH", None)
for tname, tree in TREES.items():
    if not os.path.isdir(tree):
        print(tname, "missing"); continue
    for c in cases:
        if tname == "AD" and c.endswith("get_results"):
            continue  # F-002 aborts every -fcheck run in thermo init
        with harness.slot(1):
            r = subprocess.run(["nice", "-n", "10", sys.executable, "-c", CHILD, tree, c], env=env,
                               capture_output=True, text=True, timeout=300)
        rc = r.returncode
        how = f"SIG{signal.Signals(-rc).name[3:]}" if rc < 0 else f"rc={rc}"
        out = [l for l in r.stdout.splitlines() if l.startswith(("VAL", "after", "EFTflag", "raised"))]
        err = [l for l in r.stderr.splitlines() if ".f90" in l or "Fortran runtime" in l or "Error" in l][:6]
        print(f"{tname:6s} {c:32s} {how:10s} {' | '.join(out)}")
        for l in err:
            print("        stderr:", l.strip()[:160])
