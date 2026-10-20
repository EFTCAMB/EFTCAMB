"""Track 9 (v): what does a NaN / negative / zero EFTCAMB_turn_on_time produce (no exception is raised)?

One subprocess (tree A0): Pure EFT Omega=const 0.01 (fresh CAMBparams, lmax 1000) with EFTCAMB_turn_on_time
= default, NaN, -0.1, 0, 1e-30, 2 (after today), and the GR limit (EFTflag=0). Prints TT at a few ell relative to the
default-turn-on run, plus a control "fresh_pars EFT run first, then the probe_errors reference" for state dependence.
"""

import json

from t9lib import run_snippet

CODE = r"""
def fresh_pars():
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    p.set_for_lmax(1000, lens_potential_accuracy=0)
    return p

def base_dict(**kw):
    d = {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1, 'EFTOmega0': 0.01, 'EFTwDE': 0}
    d.update(kw)
    return d

def tt(d):
    p = fresh_pars()
    p.EFTCAMB.initialize_parameters(p, d, print_header=False)
    try:
        r = camb.get_results(p)
    except Exception as e:
        return type(e).__name__ + ': ' + str(e)[:150], None
    return float(p.EFTCAMB.EFTCAMB_pert_turn_on), r.get_cmb_power_spectra(raw_cl=True, CMB_unit='muK')['total'][:, 0]

ref_on, ref = tt(base_dict())
gr_on, gr = tt({'EFTflag': 0})
out = {'default_turn_on': ref_on, 'GR_vs_default_maxrel': float(np.max(np.abs(gr[2:1000] / ref[2:1000] - 1)))}
for name, v in [('nan', float('nan')), ('negative', -0.1), ('zero', 0.0), ('tiny', 1e-30), ('after_today', 2.0)]:
    on, c = tt(base_dict(EFTCAMB_turn_on_time=v))
    if c is None:
        out[name] = on
    else:
        out[name] = {'turn_on_stored': on, 'finite': bool(np.isfinite(c[2:1000]).all()),
                     'vs_default_maxrel': float(np.max(np.abs(c[2:1000] / ref[2:1000] - 1))),
                     'vs_GR_maxrel': float(np.max(np.abs(c[2:1000] / gr[2:1000] - 1)))}
emit(**out)
"""

if __name__ == "__main__":
    r = run_snippet(CODE, threads=2, timeout=900)
    print(r["status"], json.dumps(r["result"], indent=1))
    if r["status"] != "ok":
        print("\n".join(r["tail"][-8:]))
