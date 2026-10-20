"""Track 9 (i): get_eft_functions / get_evolution / get_scale_evolution against independent routes.

Independent routes (same binary: A0 sources built with release flags + -DDEBUG, work/t9-wrapper/A0dbg):
  * the EFTCAMB debug cache files cache_*.dat written by CAMB's own source evolution (output -> derivs ->
    dump_cache_files) for every source k and time step;
  * results.hubble_parameter / conformal_time (CAMB background) for adotoa and tau;
  * the A0 reference run npz (eft_* arrays from get_eft_functions on the release tree).
Edge cases: scalar q / scalar eta, array of k vs one k at a time, unsorted eta, a before the perturbation
turn-on, a = 1, get_scale_evolution vs get_evolution(conformal_time).

Usage: python dbg_compare.py MODEL [MODEL ...]
"""

from __future__ import annotations

import json
import sys

from t9lib import AUDIT, WORK, run_snippet

FILES = {
    "FRW": "adotoa Hdot Hdotdot Hdotdotdot",
    "BDens": "grhom_t grhob_t grhoc_t grhor_t grhog_t grhov_t grhonu_tot grhonudot_tot",
    "BPres": "gpresm_t gpresdotm_t gpresdotdotm_t gpiv_t gpinu_tot gpinudot_tot gpinudotdot_tot",
    "BackgroundEFT": "EFTOmegaV EFTOmegaP EFTOmegaPP EFTOmegaPPP EFTOmegaPPPP EFTc EFTcdot EFTLambda EFTLambdadot "
    "EFTLambdadotdot",
    "SecondOrdEFT": "EFTGamma1V EFTGamma1P EFTGamma2V EFTGamma2P EFTGamma2PP EFTGamma2PPP EFTGamma3V EFTGamma3P "
    "EFTGamma3PP EFTGamma3PPP EFTGamma3PPPP EFTGamma4V EFTGamma4P EFTGamma4PP EFTGamma5V EFTGamma5P EFTGamma6V "
    "EFTGamma6P",
    "BackgroundQ": "grhoq gpresq grhodotq gpresdotq",
    "EinsteinCoeff": "EFTeomF EFTeomN EFTeomNdot EFTeomX EFTeomXdot EFTeomY EFTeomG EFTeomU EFTeomL EFTeomM EFTeomV "
    "EFTeomVdot EFTeomQ",
    "PiCoeff": "EFTpiA1 EFTpiA2 EFTpiB1 EFTpiB2 EFTpiC EFTpiD1 EFTpiD2 EFTpiE",
    "PiSolution": "pi pidot pidotdot",
    "EinsteinSol": "z sigma clxc clxb clxg clxr dgpnu dgrho dgq",
    # data order (the header of this file lists PsiDot before PhiDot, the write statement the reverse)
    "MetricMG": "Psi Phi PhiDot PsiDot mu gamma QS_mu QS_sigma QS_mu_gen QS_sigma_gen delta_DE_eff",
}
BACKGROUND_FILES = ["FRW", "BDens", "BPres", "BackgroundEFT", "SecondOrdEFT", "BackgroundQ"]
PERT_FIELDS = ["pi", "pidot", "z", "sigma", "clxc", "clxb", "clxg", "dgrho", "dgq", "Psi", "Phi", "mu", "gamma",
               "EFTeomF", "EFTeomN", "EFTeomX", "EFTeomY", "EFTpiA1", "EFTpiB1", "EFTpiC", "EFTpiE"]

CODE = r"""
import glob
FILES = {files!r}
BACKGROUND_FILES = {bgfiles!r}
PERT_FIELDS = {pfields!r}
model = PARAMS / {model!r}
d = WORK / 'dbg' / {stem!r}
d.mkdir(parents=True, exist_ok=True)
for f in glob.glob(str(d / '*')):
    os.remove(f)
wrapper = d / 'wrapper.ini'
wrapper.write_text(f'DEFAULT({{model}})\noutput_root = {{d}}/run\nEFTCAMB_write_background = F\n'
                   'l_max_scalar = 150\nget_tensor_cls = F\nget_transfer = F\ndo_lensing = F\n')
pars = camb.read_ini(str(wrapper))
results = camb.get_results(pars)
turn_on = pars.EFTCAMB.EFTCAMB_pert_turn_on

def rows(name):
    with open(d / f'run_cache_{{name}}.dat') as fh:
        for line in fh:
            t = line.strip()
            if not t or t[0] == '#' or t[0].isalpha():
                continue
            yield t

# the k column of every row, then keep only the rows of a few k values (the full dump is ~1e6 rows)
kcol = np.array([float(t.split()[2]) for t in rows('FRW')])
ks_all = np.unique(kcol)
keep_k = set([ks_all[len(ks_all) // 5], ks_all[len(ks_all) // 2], ks_all[-2]])
keep = np.array([k in keep_k for k in kcol])
def load(name):
    cols = ['a', 'tau', 'k'] + FILES[name].split()
    data = [list(map(float, t.split())) for t, m in zip(rows(name), keep) if m]
    arr = np.array(data)
    assert arr.shape[1] == len(cols), (name, arr.shape, len(cols))
    return {{c: arr[:, i] for i, c in enumerate(cols)}}

dump = {{name: load(name) for name in FILES}}
nrows_total = len(kcol)
np.savez_compressed(d / 'dump_selected.npz', **{{f'{{n}}__{{c}}': v for n in dump for c, v in dump[n].items()}})
for f in glob.glob(str(d / 'run_cache_*.dat')):
    os.remove(f)
nrows = len(dump['FRW']['a'])
for name in FILES:
    assert np.array_equal(dump[name]['a'], dump['FRW']['a']) and np.array_equal(dump[name]['k'], dump['FRW']['k'])
emit(rows_total=nrows_total, rows_kept=nrows, nk=len(ks_all), turn_on=turn_on)

def relmax(x, y):
    scale = np.max(np.abs(y))
    if scale == 0:
        return float(np.max(np.abs(x)))
    return float(np.max(np.abs(x - y)) / scale)

# ---- background: get_eft_functions at the dumped a (one k only, all its time steps) ----
k_all = dump['FRW']['k']
k0 = sorted(keep_k)[1]
sel = k_all == k0
a_d = dump['FRW']['a'][sel]
fields, ef = pars.EFTCAMB.get_eft_functions(results, a_d)
bg = {{}}
for name in BACKGROUND_FILES:
    for c in FILES[name].split():
        bg[c] = relmax(ef[c], dump[name][c][sel])
bg['tau'] = relmax(ef['tau'], dump['FRW']['tau'][sel])
emit(bg_vs_dump=bg, bg_worst=max(bg.items(), key=lambda t: t[1]))

# adotoa vs CAMB background (hubble_parameter, km/s/Mpc -> Mpc^-1 conformal)
z = 1 / a_d - 1
H = results.hubble_parameter(z) / 299792.458
emit(adotoa_vs_hubble=relmax(ef['adotoa'], a_d * H))

# vs A0 release reference (same grid as the harness)
ref = np.load(str(AUDIT_RUNS / ({stem!r} + '.npz')))
fields_g, efg = pars.EFTCAMB.get_eft_functions(results, ref['eft_a_grid'])
diffs = {{}}
for f in efg.dtype.names:
    key = 'eft_' + f
    if key in ref.files:
        diffs[f] = relmax(efg[f], ref[key])
emit(dbg_vs_A0_eftfunctions_worst=max(diffs.items(), key=lambda t: t[1]),
     dbg_vs_A0_n_nonidentical=int(sum(v != 0 for v in diffs.values())))

# ---- perturbations: get_evolution at the dumped (k, tau) for three k values ----
kpick = sorted(keep_k)
pert = {{}}
evol_single = {{}}
for kk in kpick:
    s = (k_all == kk) & (dump['FRW']['a'] > max(turn_on, 1e-5) * 1.5)
    taus = dump['FRW']['tau'][s]
    fl, ev = pars.EFTCAMB.get_evolution(results, float(kk), taus)
    evol_single[float(kk)] = (taus, ev[0])
    for f in PERT_FIELDS:
        src = next(n for n in FILES if f in FILES[n].split())
        pert.setdefault(f, []).append(relmax(ev[0][f], dump[src][f][s]))
    # the a of the evolution vs the dump
    pert.setdefault('a', []).append(relmax(ev[0]['a'], dump['FRW']['a'][s]))
emit(k_picked=[float(x) for x in kpick], pert_vs_dump={{f: max(v) for f, v in pert.items()}})

# ---- array of k at once vs one k at a time (same eta array) ----
taus = evol_single[float(kpick[1])][0]
fl, ev_multi = pars.EFTCAMB.get_evolution(results, [float(x) for x in kpick], taus)
fl, ev_one = pars.EFTCAMB.get_evolution(results, float(kpick[1]), taus)
num = [f for f in ev_one.dtype.names]
emit(multi_k_shape=list(ev_multi.shape),
     multi_vs_single_identical=all(np.array_equal(ev_multi[1][f], ev_one[0][f], equal_nan=True) for f in num))

# ---- scalar q and scalar eta ----
fl, ev_s = pars.EFTCAMB.get_evolution(results, float(kpick[1]), float(taus[-1]))
emit(scalar_shape=list(ev_s.shape), scalar_matches_last=bool(np.array_equal(ev_s[0][-1]['pi'], ev_one[0][-1]['pi'])),
     scalar_pi=float(ev_s[0][-1]['pi']), array_pi=float(ev_one[0][-1]['pi']))

# ---- unsorted eta: output order vs input order ----
t_uns = np.array([taus[-1], taus[len(taus) // 2], taus[0]])
fl, ev_u = pars.EFTCAMB.get_evolution(results, float(kpick[1]), t_uns)
emit(unsorted_input_tau=t_uns.tolist(), unsorted_output_tau=ev_u[0]['tau'].tolist(),
     unsorted_output_a=ev_u[0]['a'].tolist())

# ---- get_scale_evolution vs get_evolution(conformal_time) ----
a_req = np.array([1e-4, 0.5 * turn_on, turn_on, 2 * turn_on, 0.3, 0.7, 1.0])
fl, ev_sc = pars.EFTCAMB.get_scale_evolution(results, float(kpick[1]), a_req)
eta_req = results.conformal_time(1 / a_req - 1)
fl, ev_et = pars.EFTCAMB.get_evolution(results, float(kpick[1]), eta_req)
emit(scale_vs_eta_identical=all(np.array_equal(ev_sc[0][f], ev_et[0][f], equal_nan=True) for f in num),
     scale_a_out=ev_sc[0]['a'].tolist(), a_req=a_req.tolist(),
     pi_before_turn_on=ev_sc[0]['pi'][:2].tolist(), EFTOmegaV_before_turn_on=ev_sc[0]['EFTOmegaV'][:2].tolist(),
     EFTpiA1_before_turn_on=ev_sc[0]['EFTpiA1'][:2].tolist(), mu_after=ev_sc[0]['mu'][3:].tolist())

# ---- get_eft_functions at a = 1, scalar a, before back turn on ----
fl, e1 = pars.EFTCAMB.get_eft_functions(results, 1.0)
fl, e2 = pars.EFTCAMB.get_eft_functions(results, np.array([1.0]))
emit(eft_a1_scalar_equal=all(np.array_equal(e1[f], e2[f], equal_nan=True) for f in e1.dtype.names),
     eft_a1_tau=float(e1['tau'][0]), tau0=float(results.conformal_time(0.0)),
     eft_a1_adotoa=float(e1['adotoa'][0]), H0_conformal=float(pars.H0 / 299792.458),
     finite_a1=all(np.isfinite(e1[f]).all() for f in e1.dtype.names))
fl, e_small = pars.EFTCAMB.get_eft_functions(results, np.array([1e-12, 1e-9]))
emit(eft_tiny_a_finite={{f: bool(np.isfinite(e_small[f]).all()) for f in ('adotoa', 'EFTOmegaV', 'EFTc', 'EFTLambda',
                                                                       'tau')}},
     eft_tiny_tau=e_small['tau'].tolist())
nan_rows = {{f: int((~np.isfinite(ef[f])).sum()) for f in ef.dtype.names if (~np.isfinite(ef[f])).any()}}
emit(nonfinite_on_dump_grid=nan_rows)
"""


def main() -> None:
    models = sys.argv[1:]
    for m in models:
        stem = m.removesuffix(".ini")
        code = CODE.format(files=FILES, bgfiles=BACKGROUND_FILES, pfields=PERT_FIELDS, model=stem + ".ini",
                           stem=stem).replace("AUDIT_RUNS", f"Path({str(AUDIT / 'runs' / 'A0')!r})")
        r = run_snippet(code, tree_name=str(WORK / "A0dbg"), threads=1, timeout=3000,
                        cwd=str(WORK / "dbg"))
        out = {"model": stem, "status": r["status"], **r["result"]}
        if r["status"] != "ok":
            out["tail"] = r["tail"][-12:]
        (WORK / "dbg" / f"{stem}.result.json").write_text(json.dumps(out, indent=1))
        print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
