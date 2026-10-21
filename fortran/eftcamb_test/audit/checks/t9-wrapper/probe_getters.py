"""Probe the EFTCAMB python getters on objects without an allocated model (EFTflag = 0 / fresh CAMBparams),
and string handling (long names, non-ASCII names, values containing spaces, numpy bools, float flags)."""

import json
import sys

from t9lib import run_snippet

MODEL = "2_PEFT_Omega_power_1.ini"

PROBES = {}
for meth in ["num_params", "param_names", "param_labels", "param_values", "model_name", "get_effective_w0wa",
             "feedback", "read_parameters"]:
    PROBES[f"gr_{meth}"] = f"""
pars = camb.CAMBparams()
emit(EFTflag=pars.EFTCAMB.EFTflag)
v = pars.EFTCAMB.{meth}()
emit(value=repr(v)[:300])
"""

PROBES["set_model_name_long"] = f"""
pars = camb.read_ini(str(PARAMS / {MODEL!r}))
name = 'x' * 250 + 'END'
pars.EFTCAMB.set_model_name(name)
got = pars.EFTCAMB.model_name()
emit(len_in=len(name), len_out=len(got), truncated=(got != name))
"""
PROBES["set_model_name_utf8"] = f"""
pars = camb.read_ini(str(PARAMS / {MODEL!r}))
name = 'Omega-' + '\\u03a9' * 150
pars.EFTCAMB.set_model_name(name)
try:
    got = pars.EFTCAMB.model_name()
    emit(ok=True, got_len=len(got), equal=(got == name))
except Exception as e:
    emit(ok=False, error=repr(e)[:200])
"""
# value containing a space: helper_stripspaces removes all spaces
PROBES["value_with_space"] = """
pars = camb.CAMBparams(); pars.set_cosmology(H0=67)
pars.EFTCAMB.initialize_parameters(pars, {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': '0.0 1', 'EFTwDE': 0}, print_header=False)
emit(EFTOmega0=pars.EFTCAMB.read_parameters()['EFTOmega0'], param_values=pars.EFTCAMB.param_values())
"""
PROBES["numpy_bool"] = """
pars = camb.CAMBparams(); pars.set_cosmology(H0=67)
pars.EFTCAMB.initialize_parameters(pars, {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': 0.01, 'EFTwDE': 0, 'EFT_ghost_stability': np.True_},
                                   print_header=False)
emit(survived=True)
"""
PROBES["float_flag"] = """
pars = camb.CAMBparams(); pars.set_cosmology(H0=67)
pars.EFTCAMB.initialize_parameters(pars, {'EFTflag': 1.0, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': 0.01, 'EFTwDE': 0}, print_header=False)
emit(survived=True)
"""
PROBES["none_value"] = """
pars = camb.CAMBparams(); pars.set_cosmology(H0=67)
pars.EFTCAMB.initialize_parameters(pars, {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': None, 'EFTwDE': 0}, print_header=False)
emit(survived=True)
"""
PROBES["string_true"] = """
pars = camb.CAMBparams(); pars.set_cosmology(H0=67)
pars.EFTCAMB.initialize_parameters(pars, {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': 0.01, 'EFTwDE': 0, 'EFT_ghost_stability': 'True'},
                                   print_header=False)
emit(survived=True)
"""
PROBES["long_value_truncated"] = """
pars = camb.CAMBparams(); pars.set_cosmology(H0=67)
root = '/tmp/' + 'd' * 230 + '/out'
pars.EFTCAMB.initialize_parameters(pars, {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': 0.01, 'EFTwDE': 0, 'output_root': root}, print_header=False)
got = pars.EFTCAMB.read_parameters()['output_root']
emit(len_in=len(root), len_out=len(got), truncated=(got.decode() != root))
"""
PROBES["other_object_handle"] = f"""
import warnings
p1 = camb.read_ini(str(PARAMS / {MODEL!r}))
p2 = camb.CAMBparams(); p2.set_cosmology(H0=67)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    p1.EFTCAMB.initialize_parameters(p2, {{'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1,
                                          'EFTOmega0': 0.01, 'EFTwDE': 0}}, print_header=False)
emit(p1_flag=p1.EFTCAMB.EFTflag, p2_flag=p2.EFTCAMB.EFTflag, warnings=[str(x.message)[:120] for x in w])
"""

if __name__ == "__main__":
    sel = sys.argv[1:] or list(PROBES)
    for name in sel:
        r = run_snippet(PROBES[name], threads=1, timeout=300)
        print(f"== {name}: status={r['status']} {json.dumps(r['result'])}")
        if r["status"] != "ok":
            print("   " + " | ".join(t.strip() for t in r["tail"][-4:]))
