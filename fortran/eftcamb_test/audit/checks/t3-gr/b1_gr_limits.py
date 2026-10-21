"""t3-gr (b): GR/LCDM-equivalent points and approach-to-GR scans for each EFTCAMB model family.

python b1_gr_limits.py run [family,...]     run (cached: skips runs whose npz already exists with status ok)
python b1_gr_limits.py report [family,...]  compare with references and print/write json

For each family: GR point (default RGR), GR point with the EFT perturbation path forced on (EFTCAMB_GR_threshold=-1),
and a scan of the coupling towards zero. References (same ini, so same cosmology):
  B0   : plain CAMB v2 (EFT keys ignored)                          -> what a user compares against
  A0g  : code under audit with EFTflag=0 and IntTolBoost=2           -> isolates the EFT code path from the
         (EFTCAMB forces IntTolBoost>=2 whenever EFTflag/=0)            accuracy setting and reionization change
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from common import GR_INI, INIS, PARAMS, RUNS, WORK, dump, harness

SUB = RUNS / "b1"
P = PARAMS
SCAN = [1e-2, 1e-4, 1e-6, 1e-7, 1e-9, 0.0]

# family -> ini, GR-point overrides, scan builder(v) -> overrides, scan values, optional extra B0 reference builder
FAM: dict[str, dict] = {
    "PEFT_Omega": {"ini": P / "2_PEFT_Omega_const_1.ini", "gr": {"PureEFTmodelOmega": "0"},
                   "scan": lambda v: {"EFTOmega0": repr(v)}, "values": SCAN},
    "PEFT_Gamma1": {"ini": P / "2_PEFT_gamma1_const_1.ini", "gr": {"PureEFTmodelGamma1": "0"},
                    "scan": lambda v: {"EFTGamma10": repr(v)}, "values": [1.0, 1e-2, 1e-6, 1e-9]},
    "PEFT_wDE_LCDM": {"ini": P / "2_PEFT_wDE_LCDM_1.ini", "gr": {}, "scan": None, "values": []},
    "RPH": {"ini": P / "3_RPH_MK_const_1.ini", "gr": {"RPHmassPmodel": "0", "RPHkineticitymodel": "0"},
            "scan": lambda v: {"RPHmassP0": repr(v)}, "values": SCAN},
    "fR": {"ini": P / "4_Des_F_R_LCDM_1.ini", "gr": {"EFTB0": "0"},
           "scan": lambda v: {"EFTB0": repr(v)}, "values": [1e-2, 1e-4, 1e-6, 1e-8]},
    "mcquint": {"ini": P / "4_Des_mc5e_wCDM_1.ini", "gr": {"EFTw0": "-1.0"},
                "scan": lambda v: {"EFTw0": repr(-1.0 + v)}, "values": [1e-1, 1e-2, 1e-4, 1e-6],
                "b0_extra": lambda v: {"w": repr(-1.0 + v), "cs2_lam": "1", "dark_energy_model": "fluid"}},
    "OLgamma": {"ini": INIS / "007p5-OL-gamma.ini", "gr": {"OLOmegamodel": "0"},
                "scan": lambda v: {"OLOmega0": repr(v)}, "values": SCAN},
    "OLalpha": {"ini": INIS / "007p6-OL-alpha.ini", "gr": {"OLmassPmodel": "0", "OLkineticitymodel": "0"},
                "scan": lambda v: {"OLmass0": repr(v)}, "values": SCAN},
    "ShiftSym": {"ini": INIS / "007p7-shiftsym.ini", "gr": {"Shift_Symmetric_alphaB0": "0"},
                 "scan": lambda v: {"Shift_Symmetric_alphaB0": repr(v)}, "values": SCAN},
    "3p1": {"ini": INIS / "007p8-3p1.ini",
            "gr": {"EFTwDE": "0", "Designer3p1modelsigma1": "0", "Designer3p1modelsigma1_ODE": "0",
                   "Designer3p1modelrho1N": "0", "Designer3p1modelrho2": "0"},
            "scan": None, "values": []},
    "Horava": {"ini": P / "5_Horava_eta_1.ini", "gr": {"Horava_xi": "0", "Horava_lambda": "0", "Horava_eta": "0"},
               "scan": lambda s: {"Horava_xi": repr(-1e-5 * s), "Horava_lambda": repr(1e-5 * s),
                                  "Horava_eta": repr(2.1e-4 * s)}, "values": [1.0, 1e-2, 1e-4, 1e-6]},
    "ADE": {"ini": P / "5_ADE_1.ini", "gr": {"f_ac": "0"},
            "scan": lambda v: {"f_ac": repr(v)}, "values": [1e-2, 1e-4, 1e-6, 1e-8]},
    "Kmimic": {"ini": P / "5_Kmimic_1.ini", "gr": {"eps2_0": "0"},
               "scan": lambda v: {"eps2_0": repr(v)}, "values": [1e-2, 1e-4, 1e-6, 1e-8]},
    "Kmouflage": {"ini": P / "5_Kmouflage_1.ini", "gr": {"eps2_0": "0"},
                  "scan": lambda v: {"eps2_0": repr(-v)}, "values": [1e-2, 1e-4, 1e-6, 1e-8]},
    # Extended Galileon: B = S*q; S=0 gives EFTc = Gamma1 = Gamma2 = 0, EFTLambda = -grhov a^2 (stability needs
    # (S+0.02)*q >= 0.5, hence q=100)
    "ExtGal": {"ini": INIS / "008p7-ext-galileon.ini", "gr": {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"},
               "scan": lambda v: {"ExtendedGalileon_S": repr(v), "ExtendedGalileon_q": "100"},
               "values": [1e-2, 1e-4, 1e-6]},
    "quint_const": {"ini": P / "5_quint_quadratic.ini",
                    "gr": {"potential_model": "1", "V0": "1.0", "phidot_ini": "0"},
                    "scan": lambda v: {"p": repr(v)}, "values": [1e-1, 1e-2, 1e-4]},
    "PEFT_omk": {"ini": P / "2_PEFT_Omega_const_1.ini",
                 "gr": {"PureEFTmodelOmega": "0", "omk": "0.05", "num_redshiftwindows": "0"},
                 "scan": None, "values": []},
    "PEFT_omkm": {"ini": P / "2_PEFT_Omega_const_1.ini",
                  "gr": {"PureEFTmodelOmega": "0", "omk": "-0.05", "num_redshiftwindows": "0"},
                  "scan": None, "values": []},
    "PEFT_hmcode": {"ini": P / "2_PEFT_Omega_const_1.ini",
                    "gr": {"PureEFTmodelOmega": "0", "do_nonlinear": "2", "halofit_version": "9"},
                    "scan": None, "values": []},
    "hdsk_const": {"ini": P / "5_hdsk_TG_2.ini",
                   "gr": {"Horndeski_freefunc0_model": "1", "Horndeski_freefunc3_model": "0",
                          "Horndeski_shooting": "T", "Horndeskic00": "1.0"},
                   "scan": None, "values": []},
}


# reference cosmology groups (model inis only add EFT keys on top of these, except omnuh2 in the K-mouflage inis)
GROUPS = {"base": (GR_INI, {}), "nu0": (GR_INI, {"omnuh2": "0.0"}), "hdsk": (INIS / "t10a_GR_hdsk.ini", {}),
          "omk": (GR_INI, {"omk": "0.05", "num_redshiftwindows": "0"}),
          "omkm": (GR_INI, {"omk": "-0.05", "num_redshiftwindows": "0"}),
          "hmcode": (GR_INI, {"do_nonlinear": "2", "halofit_version": "9"})}
FAM_GROUP = {"Kmimic": "nu0", "Kmouflage": "nu0", "hdsk_const": "hdsk", "PEFT_omk": "omk", "PEFT_omkm": "omkm", "PEFT_hmcode": "hmcode"}


def ref_jobs(group: str):
    ini, ov = GROUPS[group]
    return [(f"REF_{group}__B0ref", "B0", ini, dict(ov), {}),
            (f"REF_{group}__A0g", "A0", ini, dict(ov), {"Accuracy.IntTolBoost": 2.0})]


def ref_path(fam: str, kind: str) -> Path:
    return SUB / f"REF_{FAM_GROUP.get(fam, 'base')}__{kind}.npz"


def name_of(fam: str, tag: str) -> str:
    return f"{fam}__{tag}".replace("+", "p").replace("-", "m")


def jobs_for(fam: str) -> list[tuple[str, Path, dict, dict]]:
    f = FAM[fam]
    ini = f["ini"]
    out = [
        (name_of(fam, "GR"), "A0", ini, dict(f["gr"]), {}),
        (name_of(fam, "GRforced"), "A0", ini, {**f["gr"], "EFTCAMB_GR_threshold": "-1"}, {}),
        # forced EFT path with zero EFT functions has no scalar dof: the mass-stability check then reports a spurious
        # instability, so also run it with the stability check skipped
        (name_of(fam, "GRforcedNS"), "A0", ini,
         {**f["gr"], "EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"}, {}),
    ]
    for v in f["values"]:
        out.append((name_of(fam, f"s{v:g}"), "A0", ini, f["scan"](v), {}))
        if "b0_extra" in f:
            out.append((name_of(fam, f"B0w{v:g}"), "B0", GR_INI, f["b0_extra"](v), {}))
    return out


def run_one(job) -> dict:
    name, tree, ini, ov, ps = job
    out = SUB / f"{name}.npz"
    js = out.with_suffix(".json")
    if out.exists() and js.exists() and json.loads(js.read_text()).get("status") == "ok":
        return json.loads(js.read_text())
    return harness.run_model(tree, ini, out, overrides=ov, python_set=ps, threads=1, timeout=1500)


def cmd_run(fams: list[str]) -> None:
    groups = sorted({FAM_GROUP.get(f, "base") for f in fams})
    alljobs = [j for g in groups for j in ref_jobs(g)] + [j for fam in fams for j in jobs_for(fam)]
    import os

    with ThreadPoolExecutor(int(os.environ.get("T3_JOBS", "4"))) as ex:
        for j, st in zip(alljobs, ex.map(run_one, alljobs)):
            print(f"{j[0]:40s} {st['status']:14s} {st.get('elapsed', 0):7.1f}s", flush=True)


def brief(r: dict) -> dict:
    tot = r["cl"].get("total", {})
    unl = r["cl"].get("unlensed_scalar", {})
    return {"TT": tot.get("TT"), "EE": tot.get("EE"), "TE": tot.get("TE"), "BB": tot.get("BB"),
            "uTT": unl.get("TT"), "PP": r["cl"].get("lens_potential", {}).get("PP"),
            "pk": r["worst_pk"], "s8": r.get("sigma8"), "H": r.get("bg_H"), "worst_cl": r["worst_cl"]}


def cmd_report(fams: list[str]) -> None:
    res = {}
    for fam in fams:
        res[fam] = {}
        refs = {k: ref_path(fam, k) for k in ("B0ref", "A0g")}
        for name, tree, ini, ov, ps in jobs_for(fam):
            tag = name.split("__", 1)[1]
            p = SUB / f"{name}.npz"
            st = json.loads(p.with_suffix(".json").read_text()) if p.with_suffix(".json").exists() else {}
            if st.get("status") != "ok" or not p.exists():
                res[fam][tag] = {"status": st.get("status", "missing"), "log_tail": st.get("log_tail", [])[-3:]}
                print(f"{fam:14s} {tag:12s} {st.get('status', 'missing')}")
                continue
            row = {"status": "ok"}
            for rk, rp in refs.items():
                if rp == p or not rp.exists():
                    continue
                row[rk] = brief(harness.compare(p, rp))
            if "B0w" in tag:
                pass
            elif tag.startswith("s") and "b0_extra" in FAM[fam]:
                bw = SUB / f"{name_of(fam, 'B0w' + tag[1:])}.npz"
                if bw.exists():
                    row["B0w"] = brief(harness.compare(p, bw))
            d = harness.load(p)
            row["EFTflag"] = d["meta"].get("EFTflag")
            log = p.with_suffix(".log").read_text(errors="replace")
            rgr = [ln.split(":")[-1].strip() for ln in log.splitlines() if "Return to GR time" in ln]
            row["RGR_time"] = float(rgr[-1]) if rgr else None
            res[fam][tag] = row
            line = " ".join(f"{rk}:cl={row[rk]['worst_cl']:.1e},pk={row[rk]['pk']:.1e},H={row[rk]['H'] or 0:.1e}"
                            for rk in ("B0ref", "A0g", "B0w") if rk in row)
            print(f"{fam:14s} {tag:12s} RGR={row['RGR_time']} {line}")
    dump(res, WORK / f"b1_report_{'_'.join(fams) if len(fams) < 4 else 'all'}.json")


if __name__ == "__main__":
    cmd = sys.argv[1]
    fams = sys.argv[2].split(",") if len(sys.argv) > 2 else list(FAM)
    {"run": cmd_run, "report": cmd_report}[cmd](fams)
