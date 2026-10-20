"""t8b check: background-only stability status vs H0 (cp.H0 set directly, physical densities fixed) for several
models, to see whether EFT models become 'unstable' for H0 >= 80 generically. Usage: python h0_scan_models.py TREE"""

import sys


def main(args):
    import camb

    from t8b_common import read_model

    out = {}
    for ini in args["inis"]:
        row = {}
        for H0 in args["H0s"]:
            cp = read_model(camb, ini)
            cp.H0 = H0
            try:
                camb.get_background(cp, no_thermo=True)
                row[str(H0)] = "ok"
            except Exception as e:  # noqa: BLE001
                row[str(H0)] = repr(e)[-30:-2]
        out[ini] = row
    return out


if __name__ == "__main__":
    import json

    from t8b_common import HERE, run_in_tree

    inis = ["1_EFT_GR.ini", "2_PEFT_wDE_LCDM_1.ini", "2_PEFT_Omega_power_2.ini", "4_Des_F_R_LCDM_1.ini",
            "5_quint_inverse_power.ini", "5_quint_quadratic.ini"]
    r = run_in_tree(sys.argv[1], HERE / "h0_scan_models.py", {"inis": inis, "H0s": [70, 75, 78, 80, 90]})
    for k, v in r.items():
        print(f"{k:28s} {v}")
