"""Quick probes of wrapper behaviours (each in its own subprocess): pickle/copy, write_ini, set_params(cp=...),
read_parameters conversion of a Fortran D exponent."""

import json
import sys

from t9lib import run_snippet

MODEL = "2_PEFT_Omega_power_1.ini"

PROBES = {
    "copy": f"""
pars = camb.read_ini(str(PARAMS / {MODEL!r}))
cp2 = pars.copy()
emit(flag=pars.EFTCAMB.EFTflag, flag_copy=cp2.EFTCAMB.EFTflag)
r1 = camb.get_results(pars); c1 = r1.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total']
r2 = camb.get_results(cp2); c2 = r2.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total']
emit(copy_identical=bool(np.array_equal(c1, c2)), copy_maxrel=float(np.max(np.abs(c1[2:2500,0]/c2[2:2500,0]-1))))
""",
    "pickle": f"""
import pickle
pars = camb.read_ini(str(PARAMS / {MODEL!r}))
s = pickle.dumps(pars)
p2 = pickle.loads(s)
emit(flag=pars.EFTCAMB.EFTflag, flag_unpickled=p2.EFTCAMB.EFTflag, PureEFTmodel=p2.EFTCAMB.PureEFTmodel)
try:
    emit(model_name_unpickled=p2.EFTCAMB.model_name())
except Exception as e:
    emit(model_name_error=repr(e))
r = camb.get_results(p2)
emit(got_results=True)
""",
    "write_ini_validate": f"""
pars = camb.read_ini(str(PARAMS / {MODEL!r}))
try:
    pars.write_ini(str(WORK / 'probe_written.ini'), validate=True)
    emit(write_validate='no error')
except Exception as e:
    emit(write_validate=repr(e)[:400])
pars.write_ini(str(WORK / 'probe_written_novalidate.ini'), validate=False)
txt = (WORK / 'probe_written_novalidate.ini').read_text()
emit(eft_keys_written=[l for l in txt.splitlines() if 'EFT' in l or 'Omega' in l])
p3 = camb.read_ini(str(WORK / 'probe_written_novalidate.ini'))
emit(reread_EFTflag=p3.EFTCAMB.EFTflag)
""",
    "set_params_cp": f"""
pars = camb.read_ini(str(PARAMS / {MODEL!r}))
before = pars.EFTCAMB.EFTflag
camb.set_params(cp=pars, H0=pars.H0)
emit(EFTflag_before=before, EFTflag_after=pars.EFTCAMB.EFTflag, model_after=pars.EFTCAMB.model_name())
""",
    "dexp_roundtrip": """
pars = camb.read_ini(str(PARAMS / '4_Des_F_R_LCDM_2.ini'))
d = pars.EFTCAMB.read_parameters()
emit(EFTB0=repr(d['EFTB0']), output_root=repr(d['output_root']))
cp2 = camb.read_ini(str(PARAMS / '4_Des_F_R_LCDM_2.ini'))
cp2.EFTCAMB.initialize_parameters(cp2, d, print_header=False)
emit(reinit='survived', EFTB0_after=repr(cp2.EFTCAMB.read_parameters()['EFTB0']))
""",
}

if __name__ == "__main__":
    sel = sys.argv[1:] or list(PROBES)
    for name in sel:
        r = run_snippet(PROBES[name])
        print(f"== {name}: status={r['status']}")
        print(json.dumps(r["result"], indent=1))
        if r["status"] != "ok":
            print("\n".join(r["tail"][-8:]))
