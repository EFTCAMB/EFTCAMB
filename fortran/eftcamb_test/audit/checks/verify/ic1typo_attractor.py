"""Verifier (ic1-pidot-eftpid-typo): direct check of the pi' initial condition, independent of the IC code.

For 5_Kmimic_1 we run EFTCAMB get_evolution (the real perturbation integrator) at several k on a dense time grid just
after the EFT turn-on (a_on = EFTCAMB_turn_on_time) on A0 (typo) and FIX (D2 restored). From the cached pi-equation
coefficients C, D1, D2, E (written by derivs, not by the IC routine) we rebuild in Python
    Den(a) = C + k^2 D1 + k^4 D2,   pi_att = -h0 E / Den
and the pi' the EFT_IC_type=1 formula is meant to set (with E' := 0, which is what the IC code does when any Gamma != 0,
true for K-mouflage):
    pidot_expected = h0 E Den'/Den^2,  Den' = a adotoa dDen/da (polynomial fit in ln a after a_on)
and the value the typo formula would give, with Den' -> Den' + a adotoa k^4 (dD1/da - dD2/da).
Then compare with the pidot actually integrated at a = a_on (1 + 1e-7) in each build, and measure the post-turn-on
transient |pi - pi_att| / |pi_att| up to a = 1e-6.

Usage: python ic1typo_attractor.py [turn_on_time]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

FIX = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "attractor"
KS = [1.0, 3.0, 10.0]
FIELDS = ["a", "tau", "k", "adotoa", "pi", "pidot", "EFTpiA1", "EFTpiA2", "EFTpiB1", "EFTpiB2",
          "EFTpiC", "EFTpiD1", "EFTpiD2", "EFTpiE", "EFTGamma1V", "EFTGamma2V", "EFTGamma3V", "EFTGamma4V"]


def worker(cfg):
    sys.path.insert(0, cfg["tree"])
    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(cfg["tree"]).resolve())
    cp = camb.read_ini(cfg["ini"], no_validate=True)
    res = camb.get_results(cp)  # full run (with WantCls=False the time table is not set up)
    a_on = cfg["a_on"]
    a = np.concatenate([a_on * (1 + np.array([1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 2e-2, 3e-2, 5e-2])),
                        np.logspace(np.log10(a_on * 1.1), -6, 200)])
    eta = res.conformal_time(1 / a - 1)
    fields, arr = cp.EFTCAMB.get_evolution(res, cfg["k"], eta)
    out = {f: np.ascontiguousarray(arr[f]) for f in FIELDS}
    out["a_req"] = a
    out["h0_Mpc"] = np.array(cp.H0 / 299792.458)
    np.savez(cfg["out"], **out)


def run_tree(tag, tree, a_on):
    OUT.mkdir(parents=True, exist_ok=True)
    over = {"EFTCAMB_turn_on_time": repr(a_on)}
    wrapper = harness._wrapper_ini(harness.PARAMS_DIR / "5_Kmimic_1.ini", OUT, f"{tag}_on{a_on:g}", over)
    out = OUT / f"{tag}_on{a_on:g}.npz"
    cfg = {"tree": str(tree), "ini": str(wrapper), "out": str(out), "k": KS, "a_on": a_on}
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "_worker", json.dumps(cfg)], cwd=OUT, env=env,
                           capture_output=True, text=True, timeout=900, preexec_fn=harness._nice)
    if p.returncode != 0:
        raise RuntimeError(f"{tag}: rc={p.returncode}\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    with np.load(out) as f:
        return {k: f[k] for k in f.files}


def analyse(d, a_on):
    h0 = float(d["h0_Mpc"])
    res = {}
    for ik, k in enumerate(KS):
        g = lambda f: d[f][ik, :]  # noqa: E731  get_evolution returns (nq, neta)
        a = g("a")
        assert np.allclose(g("k"), k) and np.allclose(a, d["a_req"], rtol=1e-6), (g("k")[:3], a[:3], d["a_req"][:3])
        C, D1, D2, E, H = g("EFTpiC"), g("EFTpiD1"), g("EFTpiD2"), g("EFTpiE"), g("adotoa")
        Den = C + k**2 * D1 + k**4 * D2
        pi_att = -h0 * E / Den
        # derivatives at a_on from the first 12 points (a_on*(1+1e-7) .. a_on*1.05): quadratic fit in x = ln(a/a_on)
        sel = slice(0, 12)
        x = np.log(a[sel] / a_on)
        dDen_dx = np.polyfit(x, Den[sel], 3)[-2]
        dD1_dx = np.polyfit(x, D1[sel], 3)[-2]
        dD2_dx = np.polyfit(x, D2[sel], 3)[-2]
        i0 = 0
        # d/dtau = a adotoa d/da = adotoa d/dx
        Denp = H[i0] * dDen_dx
        typo_extra = H[i0] * k**4 * (dD1_dx - dD2_dx)
        pd_exp = h0 * E[i0] * Denp / Den[i0] ** 2
        pd_typo = h0 * E[i0] * (Denp + typo_extra) / Den[i0] ** 2
        late = a > 2 * a_on
        dev = np.abs(g("pi")[late] - pi_att[late]) / np.abs(pi_att[late])
        res[k] = {
            "a0": float(a[i0]), "pi0": float(g("pi")[i0]), "pi_att0": float(pi_att[i0]),
            "pidot0": float(g("pidot")[i0]), "pidot_expected_fixedformula": float(pd_exp),
            "pidot_expected_typoformula": float(pd_typo),
            "k2_D1prime_over_Cprime_plus": float(k**2 * dD1_dx / dDen_dx),
            "D2_at_on": float(D2[i0]), "Gamma1": float(g("EFTGamma1V")[i0]),
            "max_rel_dev_pi_vs_att_after_2aon": float(np.max(dev)),
            "median_rel_dev": float(np.median(dev)),
        }
    return res


def main():
    a_on = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-8
    out = {}
    for tag, tree in (("A0", harness.TREES["A0"]), ("FIX", FIX)):
        d = run_tree(tag, tree, a_on)
        out[tag] = analyse(d, a_on)
    print(f"5_Kmimic_1, turn-on a_on = {a_on:g}")
    for k in KS:
        a, f = out["A0"][k], out["FIX"][k]
        print(f"k = {k:g}/Mpc  (D2 = {a['D2_at_on']:.3e}, Gamma1 = {a['Gamma1']:.3e})")
        print(f"   pi at a_on(1+1e-7):  A0 {a['pi0']:+.6e}  FIX {f['pi0']:+.6e}  pi_att {a['pi_att0']:+.6e}")
        print(f"   pidot integrated:    A0 {a['pidot0']:+.6e}  FIX {f['pidot0']:+.6e}")
        print(f"   pidot from cached coefficients (independent): correct formula {a['pidot_expected_fixedformula']:+.6e}"
              f"   typo formula {a['pidot_expected_typoformula']:+.6e}")
        print(f"   |pi-pi_att|/|pi_att| for a in (2a_on, 1e-6]: max A0 {a['max_rel_dev_pi_vs_att_after_2aon']:.3e} "
              f"FIX {f['max_rel_dev_pi_vs_att_after_2aon']:.3e};  median A0 {a['median_rel_dev']:.3e} "
              f"FIX {f['median_rel_dev']:.3e}")
    (OUT / f"summary_on{a_on:g}.json").write_text(json.dumps({t: {str(k): v for k, v in r.items()}
                                                             for t, r in out.items()}, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        worker(json.loads(sys.argv[2]))
    else:
        main()
