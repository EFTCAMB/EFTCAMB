"""Run named sets of (tree, ini, overrides, python_set) for track 10a, skipping outputs that already exist.

Usage: python run_sets.py SETNAME [--jobs 4] [--models a,b,c] [--force]

Sets:
  grrefs      GR references (EFTflag=0) at IntTolBoost 1 and 2 for base, omnuh2=0 and hdsk bases, A0 and B1
  ha          all models + GR refs at accuracy_boost=2, IntTolBoost=2 (python_set), A0 and B1
  itb2_b1     all models in B1 with IntTolBoost=2 (python_set) at default accuracy
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = Path(harness.AUDIT_ROOT) / "work" / "t10a-regression"
INIS = HERE.parents[1] / "inis"
GR_INIS = {
    "1_EFT_GR": harness.PARAMS_DIR / "1_EFT_GR.ini",
    "t10a_GR_nu0": INIS / "t10a_GR_nu0.ini",
    "t10a_GR_hdsk": INIS / "t10a_GR_hdsk.ini",
}
WCDM_INI = INIS / "t10a_GR_wCDM08.ini"
CONV_MODELS = ["2_PEFT_Omega_const_1", "2_PEFT_gamma3_const_1", "2_PEFT_gamma4_const_1", "2_PEFT_wDE_wCDM_1",
               "4_Des_F_R_LCDM_1", "5_Kmouflage_2", "5_quint_cubic", "3_RPH_MK_const_1"]


def jobs_for(setname: str, models: list[str] | None) -> list[tuple]:
    out = []
    allm = {p.stem: p for p in harness.model_inis()}
    if models:
        allm = {k: v for k, v in allm.items() if k in models}
    if setname == "grrefs":
        for tree in ("A0", "B1"):
            for itb in (1, 2):
                for name, ini in GR_INIS.items():
                    ps = {"Accuracy.IntTolBoost": float(itb)} if itb != 1 else {}
                    out.append((tree, ini, WORK / f"runs_{tree}_itb{itb}" / f"{name}.npz", {}, ps))
    elif setname == "ha":
        ov = {"accuracy_boost": "2"}
        ps = {"Accuracy.IntTolBoost": 2.0}
        for tree in ("A0", "B1"):
            for name, ini in {**GR_INIS, **allm}.items():
                out.append((tree, ini, WORK / f"runs_{tree}_ha" / f"{name}.npz", ov, ps))
    elif setname == "itb2_b1":
        ps = {"Accuracy.IntTolBoost": 2.0}
        for name, ini in allm.items():
            out.append(("B1", ini, WORK / "runs_B1_itb2" / f"{name}.npz", {}, ps))
    elif setname == "nu0":
        # massless neutrinos only (omnuh2 = 0), IntTolBoost = 2 in both trees (A0 forces it for EFT models)
        ov = {"omnuh2": "0.0"}
        ps = {"Accuracy.IntTolBoost": 2.0}
        for tree in ("A0", "B1"):
            for name, ini in {"1_EFT_GR": GR_INIS["1_EFT_GR"], "t10a_GR_hdsk": GR_INIS["t10a_GR_hdsk"], **allm}.items():
                out.append((tree, ini, WORK / f"runs_{tree}_nu0" / f"{name}.npz", ov, ps))
    elif setname == "wcdm2":
        refs = {n: INIS / f"{n}.ini" for n in ("t10a_GR_wCDM09", "t10a_GR_CPL_09_01", "t10a_GR_CPL_07_03")}
        for tree, ps, tag in (("A0", {"Accuracy.IntTolBoost": 2.0}, "itb2"), ("B1", {}, "itb1")):
            for name, ini in refs.items():
                out.append((tree, ini, WORK / f"runs_{tree}_{tag}" / f"{name}.npz", {}, ps))
    elif setname == "wtab":
        refs = {x.stem: x for x in sorted(INIS.glob("t10a_GR_w_*.ini"))}
        for tree, ps, tag in (("A0", {"Accuracy.IntTolBoost": 2.0}, "itb2"), ("B1", {}, "itb1")):
            for name, ini in refs.items():
                out.append((tree, ini, WORK / f"runs_{tree}_{tag}" / f"{name}.npz", {}, ps))
    elif setname == "ab2i1":
        # B1 at accuracy_boost=2 with its default IntTolBoost=1 (B1 dverk fails with ab2 + IntTolBoost=2)
        sel = {"1_EFT_GR": GR_INIS["1_EFT_GR"], "t10a_GR_nu0": GR_INIS["t10a_GR_nu0"], "t10a_GR_wCDM08": WCDM_INI}
        sel.update({m: allm[m] for m in CONV_MODELS if m in allm})
        for name, ini in sel.items():
            out.append(("B1", ini, WORK / "runs_B1_ab2i1" / f"{name}.npz", {"accuracy_boost": "2"}, {}))
    elif setname == "weff":
        # GR references with the GR-equivalent dark energy of MG expansion histories (make_wtab_refs.py --eff)
        refs = {x.stem: x for x in sorted(INIS.glob("t10a_GR_weff_*.ini"))}
        for tree, ps, tag in (("A0", {"Accuracy.IntTolBoost": 2.0}, "itb2"), ("B1", {}, "itb1")):
            for name, ini in refs.items():
                out.append((tree, ini, WORK / f"runs_{tree}_{tag}" / f"{name}.npz", {}, ps))
    elif setname == "lab2sel":
        # l_accuracy_boost = 2 (photon/neutrino hierarchy lmax) in both trees for selected models
        sel = ["2_PEFT_gamma3_const_1", "3_RPH_MK_const_1", "4_Des_F_R_LCDM_1", "2_PEFT_Omega_lin_1"]
        for name in sel:
            for tree in ("B1", "A0"):
                ps = {"Accuracy.IntTolBoost": 2.0} if tree == "A0" else {}
                out.append((tree, allm[name], WORK / f"runs_{tree}_lab2" / f"{name}.npz", {"l_accuracy_boost": "2"}, ps))
    elif setname == "nu0ab2":
        # massless neutrinos, accuracy_boost = 2, IntTolBoost = 2 in both trees (B1 dverk fails at ab2 with massive nu)
        ov = {"omnuh2": "0.0", "accuracy_boost": "2"}
        ps = {"Accuracy.IntTolBoost": 2.0}
        # (1_EFT_GR with omnuh2=0 at ab2 is identical to runs_*_ab2/t10a_GR_nu0: linked, not rerun)
        sel = ["2_PEFT_Omega_const_1", "3_RPH_MK_const_1", "2_PEFT_gamma3_const_1", "4_Des_F_R_LCDM_1",
               "2_PEFT_gamma4_const_1", "2_PEFT_gamma5_const_1", "2_PEFT_Omega_lin_1", "5_Horava_eta_4",
               "3_RPH_MB_OmegaDE_1"]
        for name in sel:
            for tree in ("B1", "A0"):
                out.append((tree, allm[name], WORK / f"runs_{tree}_nu0ab2" / f"{name}.npz", ov, ps))
    elif setname == "conv_full":
        # converged comparison: massless nu, accuracy_boost=2, l_accuracy_boost=2, IntTolBoost=2 in both trees
        ov = {"omnuh2": "0.0", "accuracy_boost": "2", "l_accuracy_boost": "2"}
        ps = {"Accuracy.IntTolBoost": 2.0}
        sel = {"1_EFT_GR": GR_INIS["1_EFT_GR"]}
        sel.update({m: allm[m] for m in ("2_PEFT_Omega_const_1", "2_PEFT_gamma3_const_1", "2_PEFT_gamma4_const_1",
                                         "3_RPH_MK_const_1", "5_Horava_eta_4", "5_scg", "4_Des_F_R_LCDM_1")})
        for name, ini in sel.items():
            for tree in ("B1", "A0"):
                out.append((tree, ini, WORK / f"runs_{tree}_conv" / f"{name}.npz", ov, ps))
    elif setname == "hdsk":
        # Horndeski (hdsk base): massless nu at default accuracy and at converged settings, both trees (IntTolBoost=2)
        ps = {"Accuracy.IntTolBoost": 2.0}
        sel = {"t10a_GR_hdsk": GR_INIS["t10a_GR_hdsk"], "5_hdsk_TG_1": allm["5_hdsk_TG_1"]}
        for tag, ov in (("nu0", {"omnuh2": "0.0"}),
                        ("conv", {"omnuh2": "0.0", "accuracy_boost": "2", "l_accuracy_boost": "2"})):
            for name, ini in sel.items():
                for tree in ("B1", "A0"):
                    out.append((tree, ini, WORK / f"runs_{tree}_{tag}" / f"{name}.npz", ov, ps))
    elif setname == "kmimic_ab2":
        # massless-neutrino K-mimic / K-mouflage models at accuracy_boost=2, IntTolBoost=2 (B1 works without massive nu)
        for name in ("5_Kmimic_1", "5_Kmouflage_1"):
            for tree in ("B1", "A0"):
                out.append((tree, allm[name], WORK / f"runs_{tree}_ab2" / f"{name}.npz", {"accuracy_boost": "2"},
                            {"Accuracy.IntTolBoost": 2.0}))
    elif setname == "b1fix":
        # B1 copy whose ThermalNuBackground%pidot/pidotdot/pidotdotdot are replaced by the A0 ones (patch_b1_nufix.py)
        tree = str(WORK / "B1nufix")
        prio = ["1_EFT_GR", "5_scg", "2_PEFT_Omega_const_1", "3_RPH_MK_const_1", "4_Des_F_R_LCDM_1",
                "2_PEFT_gamma3_const_1", "2_PEFT_gamma4_const_1", "5_Horava_eta_4", "5_quint_cubic", "5_ADE_1",
                "t10a_GR_hdsk", "5_hdsk_TG_1", "4_Des_mc5e_turning_point_2", "2_PEFT_Omega_lin_1"]
        names = prio + sorted(m for m in allm if m not in prio)
        refs = {"1_EFT_GR": GR_INIS["1_EFT_GR"], "t10a_GR_hdsk": GR_INIS["t10a_GR_hdsk"]}
        for name in names:
            ini = refs.get(name) or allm.get(name)
            if ini is None or name.startswith(("5_Kmouflage", "5_Kmimic")):
                continue
            out.append((tree, ini, WORK / "runs_B1fix_itb1" / f"{name}.npz", {}, {}))
    elif setname == "nu0sel":
        # selected massive-nu models rerun with omnuh2 = 0 in both trees (IntTolBoost = 2 in both)
        ov = {"omnuh2": "0.0"}
        ps = {"Accuracy.IntTolBoost": 2.0}
        sel = ["2_PEFT_Omega_const_1", "2_PEFT_Omega_lin_1", "2_PEFT_gamma3_const_1", "2_PEFT_gamma4_const_1",
               "2_PEFT_gamma5_const_1", "3_RPH_MK_const_1", "3_RPH_MB_OmegaDE_1", "4_Des_F_R_LCDM_1",
               "4_Des_F_R_wCDM_1", "4_Des_mc5e_turning_point_2", "5_ADE_1", "5_Horava_eta_4", "5_quint_cubic",
               "5_quint_axion", "5_scg", "2_PEFT_wDE_wCDM_1"]
        for tree in ("A0", "B1"):
            for name, ini in {"1_EFT_GR": GR_INIS["1_EFT_GR"], "t10a_GR_wCDM08": WCDM_INI,
                              **{m: allm[m] for m in sel}}.items():
                out.append((tree, ini, WORK / f"runs_{tree}_nu0" / f"{name}.npz", ov, ps))
    elif setname == "noreion":
        # no reionization, massless neutrinos, IntTolBoost = 2 in both trees
        ov = {"omnuh2": "0.0", "reionization": "F"}
        ps = {"Accuracy.IntTolBoost": 2.0}
        for tree in ("A0", "B1"):
            for name, ini in {"1_EFT_GR": GR_INIS["1_EFT_GR"], "t10a_GR_hdsk": GR_INIS["t10a_GR_hdsk"], **allm}.items():
                out.append((tree, ini, WORK / f"runs_{tree}_noreion" / f"{name}.npz", ov, ps))
    elif setname in ("conv_itb4", "conv_ab2", "wcdm"):
        sel = {"1_EFT_GR": GR_INIS["1_EFT_GR"], "t10a_GR_nu0": GR_INIS["t10a_GR_nu0"], "t10a_GR_wCDM08": WCDM_INI}
        if setname != "wcdm":
            sel.update({m: allm[m] for m in CONV_MODELS if m in allm})
        if setname == "conv_itb4":
            cfgs = [("A0", {}, {"Accuracy.IntTolBoost": 4.0}, "itb4"), ("B1", {}, {"Accuracy.IntTolBoost": 4.0}, "itb4")]
        elif setname == "conv_ab2":
            cfgs = [(t, {"accuracy_boost": "2"}, {"Accuracy.IntTolBoost": 2.0}, "ab2") for t in ("A0", "B1")]
        else:
            cfgs = [("A0", {}, {"Accuracy.IntTolBoost": 2.0}, "itb2"), ("B1", {}, {}, "itb1"),
                    ("B1", {}, {"Accuracy.IntTolBoost": 2.0}, "itb2"), ("A0", {}, {}, "itb1")]
        for tree, ov, ps, tag in cfgs:
            for name, ini in sel.items():
                out.append((tree, ini, WORK / f"runs_{tree}_{tag}" / f"{name}.npz", ov, ps))
    else:
        raise SystemExit(f"unknown set {setname}")
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("setname")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--models", default="")
    p.add_argument("--force", action="store_true")
    p.add_argument("--trees", default="A0,B1")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--budget", type=float, default=0.0, help="seconds after which no new run is started (0 = none)")
    p.add_argument("--max-runs", type=int, default=0)
    a = p.parse_args()
    models = a.models.split(",") if a.models else None
    todo = [j for j in jobs_for(a.setname, models) if j[0] in a.trees.split(",")]
    if not a.force:
        todo = [j for j in todo if not j[2].exists()]
    if a.max_runs:
        todo = todo[: a.max_runs]
    print(f"{len(todo)} runs", flush=True)
    import time

    t0 = time.time()

    def one(j):
        tree, ini, out, ov, ps = j
        if a.budget and time.time() - t0 > a.budget:
            print(f"{tree} {out.parent.name}/{out.stem:34s} skipped (budget)", flush=True)
            return None
        r = harness.run_model(tree, ini, out, overrides=ov, python_set=ps, threads=a.threads, timeout=3000)
        print(f"{tree} {out.parent.name}/{out.stem:34s} {r['status']:16s} {r['elapsed']:7.1f}s", flush=True)
        return r

    with ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(one, todo))


if __name__ == "__main__":
    main()
