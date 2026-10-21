"""t8b check: camb/model.py set_cosmology EFTCAMB MOD restricting the H0 bracket of set_H0_for_theta to
H0 > 1.1*100*sqrt(omega_b + omega_nu + omega_nu_sterile + omega_c), applied also to GR.
Closed models with a low H0 solution become unsolvable. Usage: python theta_bracket.py TREE..."""

import sys


def main(args):
    import camb

    out = {}
    for omk in (0.0, -0.05, -0.1, -0.15):
        for meff in (0.0,):
            try:
                pars = camb.set_params(cosmomc_theta=0.0104109, ombh2=0.0224, omch2=0.12, omk=omk, mnu=0.06,
                                       As=2.1e-9, ns=0.965)
                out[f"omk={omk}"] = {"H0": pars.H0}
            except Exception as e:  # noqa: BLE001
                out[f"omk={omk}"] = {"error": repr(e)[:160]}
    # sterile neutrino: is omnuh2_sterile double counted in min_H0? (self.omnuh2 already includes it)
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.0, ombh2=0.0224, omch2=0.12, mnu=0.06, nnu=3.6, meffsterile=0.5)
    out["sterile_omnuh2_total"] = pars.omnuh2
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    for t in sys.argv[1:]:
        print(t, json.dumps(run_in_tree(t, HERE / "theta_bracket.py", {}), indent=1))
