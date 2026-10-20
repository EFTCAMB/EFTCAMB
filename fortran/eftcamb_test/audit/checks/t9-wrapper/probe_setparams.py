"""Track 9 (iii)+(vii): camb.set_params interplay with EFTCAMB.

Probes (tree A0, each in its own subprocess):
  cp_reset     : set_params(cp=<EFT params from read_ini>, H0=<same H0>) -> is the EFT model kept?
  de_leak      : set_params(**eft_dict) vs set_params(**eft_dict, w=-0.8, wa=0.1): does the CAMB DarkEnergy w
                 change an EFT model (double counting / leak)?
  collisions   : EFT keys read by any test model (from the roundtrip json files) vs every keyword set_params
                 accepts (setter signatures + CAMBparams attributes); and which CAMB keys EFTCAMB reads.
  camb_keys_in_eft : set_params(**eft_dict, lmax=1200, ...) -> which non-EFT keys end up in read_parameters().

Usage: python probe_setparams.py [probe ...]
"""

from __future__ import annotations

import json
import sys

from t9lib import WORK, run_snippet

MODEL = "2_PEFT_wDE_wCDM_1.ini"

COMMON = r"""
def eft_dict(ini):
    p = camb.read_ini(str(PARAMS / ini))
    return p, {k: (v.decode() if isinstance(v, bytes) else v) for k, v in p.EFTCAMB.read_parameters().items()}

def cls(p):
    r = camb.get_results(p)
    return r.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:1200, 0]

def base(**kw):
    d = dict(H0=67.5, ombh2=0.022, omch2=0.12, As=2.1e-9, ns=0.965, lmax=1200, lens_potential_accuracy=0)
    d.update(kw)
    return d
"""

PROBES = {
    "cp_reset": f"""
p, d = eft_dict({MODEL!r})
p.set_for_lmax(1200, lens_potential_accuracy=0); p.WantTensors = False
before = (int(p.EFTCAMB.EFTflag), p.EFTCAMB.model_name(), p.EFTCAMB.param_values())
c0 = cls(p)
camb.set_params(cp=p, H0=p.H0)
after = (int(p.EFTCAMB.EFTflag), p.EFTCAMB.model_name(), p.EFTCAMB.param_values())
c1 = cls(p)
emit(before=before, after=after, TT_maxrel_change=float(np.max(np.abs(c1 / c0 - 1))))
""",
    "de_leak": f"""
p0, d = eft_dict({MODEL!r})
d.pop('output_root', None)
pa = camb.set_params(**base(**d))
pb = camb.set_params(**base(w=-0.8, wa=0.1, **d))
pc = camb.set_params(**base(dark_energy_model='ppf', w=-0.8, wa=0.1, **d))
emit(model=pa.EFTCAMB.model_name(), w_a=pa.DarkEnergy.w, w_b=pb.DarkEnergy.w, de_class_c=type(pc.DarkEnergy).__name__)
ca, cb, cc = cls(pa), cls(pb), cls(pc)
ra = camb.get_results(pa); rb = camb.get_results(pb)
emit(TT_maxrel_w_fluid=float(np.max(np.abs(cb / ca - 1))), TT_maxrel_w_ppf=float(np.max(np.abs(cc / ca - 1))),
     H_z1_rel=float(rb.hubble_parameter(1.0) / ra.hubble_parameter(1.0) - 1))
# GR control: the same w in plain CAMB must change the spectra
g1 = camb.set_params(**base()); g2 = camb.set_params(**base(w=-0.8, wa=0.1))
emit(GR_TT_maxrel_w=float(np.max(np.abs(cls(g2) / cls(g1) - 1))))
""",
    "collisions": r"""
import glob, inspect
from camb.camb import _get_setter_spec
cp = camb.CAMBparams()
setters = [cp.set_accuracy, cp.set_classes, cp.Recomb.set_params, cp.DarkEnergy.set_params, cp.Reion.set_extra_params,
           cp.set_cosmology, cp.set_matter_power, cp.set_for_lmax, cp.InitPower.set_params,
           cp.NonLinearModel.set_params]
camb_keys = set()
for s in setters:
    camb_keys |= set(_get_setter_spec(s))
attr_keys = {n for n in dir(cp) if not n.startswith('_')}
eft_keys = set()
for f in glob.glob(str(WORK / 'roundtrip' / '*.json')):
    eft_keys |= set(json.load(open(f)).get('types', {}))
emit(n_eft_keys=len(eft_keys), n_camb_setter_keys=len(camb_keys),
     eft_vs_setters=sorted(eft_keys & camb_keys), eft_vs_attributes=sorted(eft_keys & attr_keys),
     eft_keys_case_insensitive_vs_camb=sorted(k for k in eft_keys if k.lower() in {c.lower() for c in camb_keys | attr_keys}),
     non_eft_prefixed=sorted(k for k in eft_keys if not any(t in k for t in ('EFT', 'RPH', 'Horndeski', 'Horava',
         'Kmimic', 'Kmouflage', 'Designer', 'Shift', 'Scaling', 'Beyond', 'Quint', 'quint', 'ADE', 'Galileon')))[:80])
""",
    "camb_keys_in_eft": f"""
p0, d = eft_dict({MODEL!r})
d.pop('output_root', None)
p = camb.set_params(**base(**d), WantTransfer=True, redshifts=[0.0, 1.0], kmax=2.0, Alens=1.1)
rp = p.EFTCAMB.read_parameters()
emit(n_read=len(rp), non_eft_read=sorted(k for k in rp if k not in d))
""",
    # 'cs2' is both a CAMB set_params key (DarkEnergy.set_params) and the ADE sound speed key (008p2_ADE.f90:181)
    "ade_cs2": r"""
f = WORK / 'ade_cs2.ini'
f.write_text(f'DEFAULT({PARAMS}/5_ADE_1.ini)\ncs2 = 0.5\nget_tensor_cls = F\n')
p1 = camb.read_ini(str(f))
d = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in p1.EFTCAMB.read_parameters().items()}
d.pop('output_root', None)
g = camb.read_ini(str(f)); g.EFTCAMB.initialize_parameters(g, {'EFTflag': 0}, print_header=False)
p2 = camb.set_params(cp=g, **d)
c1, c2 = cls(p1), cls(p2)
emit(cs2_read=d.get('cs2'), DE_cs2_ini=float(p1.DarkEnergy.cs2), DE_cs2_setparams=float(p2.DarkEnergy.cs2),
     eft_cs2_setparams=p2.EFTCAMB.read_parameters().get('cs2'), TT_maxrel=float(np.max(np.abs(c2 / c1 - 1))))
""",
    # a long, non-scalar CAMB value is also written into the EFT ini reader
    "long_camb_value": f"""
p0, d = eft_dict({MODEL!r})
d.pop('output_root', None)
z = list(np.linspace(0.0, 5.0, 60)[::-1])
pa = camb.set_params(**base(**d), WantTransfer=True, redshifts=z, kmax=1.0)
pb = camb.set_params(**base(**d), WantTransfer=True, redshifts=z[:2], kmax=1.0)
emit(len_str=len(str(z)), same_eft=(pa.EFTCAMB.read_parameters() == pb.EFTCAMB.read_parameters()),
     model=pa.EFTCAMB.model_name())
""",
}

if __name__ == "__main__":
    sel = sys.argv[1:] or list(PROBES)
    for name in sel:
        r = run_snippet(COMMON + PROBES[name], threads=2, timeout=900)
        print(f"== {name}: status={r['status']} {json.dumps(r['result'], default=str)}", flush=True)
        if r["status"] != "ok":
            print("   " + " | ".join(t.strip() for t in r["tail"][-6:]))
