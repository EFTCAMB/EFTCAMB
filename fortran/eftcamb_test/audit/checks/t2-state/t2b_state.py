"""Track 2b: state leakage between successive EFTCAMB calculations in one Python process.

A *step* computes one model and saves a fingerprint npz (Cls, P(k), sigma8, H(z), EFT functions, derived).  A
*sequence* runs several steps in one subprocess (OMP_NUM_THREADS=1).  Every step is also run alone in a fresh process.
Each step's output in a sequence must be bit identical to the fresh-process output and to its earlier occurrence.

Step kinds (dicts):
  {"kind": "ini", "model": M, "overrides": {...}}          read_ini(DEFAULT(M.ini)+overrides) + get_results
  {"kind": "setparams", "model": M, "cosmo": {...}}         camb.set_params(**cosmo, **EFT keys of M.ini)
  {"kind": "reinit", "model": M, "cosmo": {...}}            reuse the CAMBparams of the previous step and call
                                                            cp.EFTCAMB.initialize_parameters(cp, EFT keys of M)
  {"kind": "deepcopy" | "copy" | "pickle", "model": M, ...} set_params(M), then compute from copy.deepcopy(cp) /
                                                            cp.copy() / pickle.loads(pickle.dumps(cp))
Usage:
  python t2b_state.py all            run every sequence defined in SEQUENCES and print the comparison
  python t2b_state.py seq NAME       run one sequence
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "t2-state" / "state"
TREE = os.environ.get("T2B_TREE", "A0")
A_GRID = np.logspace(-6, 0, 200)

# low-cost settings shared by all steps (keeps the sequences fast; the same settings are used in the fresh runs)
FAST = {"l_max_scalar": "800", "l_max_tensor": "400", "k_eta_max_tensor": "800", "feedback_level": "0"}
COSMO = dict(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, tau=0.09, As=2.1e-9, ns=0.96, lmax=800, r=0.1,
             WantTensors=True, feedback_level=0)
COSMO_NONU = dict(COSMO, mnu=0.0, num_massive_neutrinos=0)
NONU = {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046"}
UNSTABLE = {"EFTOmega0": "-1.0"}
BADFLAG = {"PureEFTmodel": "7"}


def ini(model, **ov):
    return {"kind": "ini", "model": model, "overrides": {**FAST, **{k: str(v) for k, v in ov.items()}}}


def sp(model, kind="setparams", cosmo=None, **eft):
    st = {"kind": kind, "model": model, "cosmo": cosmo or COSMO}
    if eft:
        st["eft"] = {k: str(v) for k, v in eft.items()}
    return st


SEQUENCES = {
    # (b1) EFTflag 0 -> nonzero -> 0
    "gr_eft_gr": [ini("1_EFT_GR"), ini("2_PEFT_Omega_const_1"), ini("1_EFT_GR")],
    # (b2) designer -> full mapping -> designer
    "des_full_des": [ini("4_Des_F_R_LCDM_1"), ini("5_Kmouflage_3"), ini("4_Des_F_R_LCDM_1")],
    "mc5e_quint_mc5e": [ini("4_Des_mc5e_CPL_1"), ini("5_quint_exponential"), ini("4_Des_mc5e_CPL_1")],
    # (b3) two full-mapping models
    "full_full": [ini("5_quint_exponential"), ini("5_ADE_1"), ini("5_Horava_eta_1"), ini("5_quint_exponential"),
                  ini("5_ADE_1")],
    "kmimic_kmouflage": [ini("5_Kmimic_1"), ini("5_Kmouflage_1"), ini("5_Kmimic_1")],
    # (b4) two parametrisations of the same model
    "peft_param": [ini("2_PEFT_Omega_const_1"), ini("2_PEFT_Omega_exp_1"), ini("2_PEFT_gamma3_lin_1"),
                   ini("2_PEFT_Omega_const_1")],
    "rph_param": [ini("3_RPH_KT_OmegaDE_1"), ini("3_RPH_MK_power_1"), ini("3_RPH_KT_OmegaDE_1")],
    "mc5e_param": [ini("4_Des_mc5e_wCDM_1"), ini("4_Des_mc5e_taylor_1"), ini("4_Des_mc5e_wCDM_1")],
    "fr_param": [ini("4_Des_F_R_wCDM_2"), ini("4_Des_F_R_LCDM_2"), ini("4_Des_F_R_wCDM_2")],
    # (b5) with / without massive neutrinos
    "nu_nonu": [ini("2_PEFT_Omega_const_1"), ini("2_PEFT_Omega_const_1", **NONU), ini("2_PEFT_Omega_const_1"),
                ini("4_Des_F_R_LCDM_1", **NONU), ini("4_Des_F_R_LCDM_1"), ini("4_Des_F_R_LCDM_1", **NONU)],
    # (b6) failing model between two good ones (unstable theory, invalid model flag)
    "fail_unstable": [ini("2_PEFT_Omega_const_1"), ini("2_PEFT_Omega_const_1", **UNSTABLE),
                      ini("2_PEFT_Omega_const_1")],
    "fail_badflag": [ini("2_PEFT_Omega_exp_1"), ini("2_PEFT_Omega_const_1", **BADFLAG), ini("2_PEFT_Omega_exp_1")],
    "fail_unstable_gr": [ini("1_EFT_GR"), ini("2_PEFT_Omega_const_1", **UNSTABLE), ini("1_EFT_GR")],
    # set_params path
    "sp_seq": [sp("2_PEFT_Omega_const_1"), sp("5_quint_exponential"), sp("4_Des_F_R_LCDM_1"),
               sp("2_PEFT_Omega_const_1"), sp("1_EFT_GR"), sp("5_quint_exponential"),
               sp("4_Des_F_R_LCDM_1", cosmo=COSMO_NONU), sp("4_Des_F_R_LCDM_1")],
    # reuse one CAMBparams object, re-initialising EFTCAMB in place
    "reinit": [sp("2_PEFT_Omega_const_1"), sp("4_Des_F_R_LCDM_1", "reinit"), sp("2_PEFT_Omega_const_1", "reinit"),
               sp("5_quint_exponential", "reinit"), sp("1_EFT_GR", "reinit"), sp("2_PEFT_Omega_exp_1", "reinit")],
    # Horndeski: SAVEd num_eq (evolve_hubble T -> F) and the never-filled V_scf table (fingerprint has eft_V_scf)
    "hdsk_seq": [ini("5_hdsk_TG_1"), ini("5_quint_exponential"), ini("5_hdsk_TG_1"), ini("5_hdsk_TG_3"),
                 ini("5_hdsk_TG_1")],
    "hdsk_hubble": [ini("5_hdsk_TG_1"), ini("5_hdsk_TG_1", Horndeski_evolve_hubble="T"), ini("5_hdsk_TG_1")],
    # same CAMBparams object used for two get_results calls (no re-initialisation)
    "samecp": [sp("4_Des_F_R_LCDM_1", "samecp"), sp("5_quint_exponential", "samecp")],
    # one CAMBdata object reused: results.calc_power_spectra(params) with different models
    "samedata": [sp("2_PEFT_Omega_const_1", "samedata"), sp("4_Des_F_R_LCDM_1", "samedata"),
                 sp("2_PEFT_Omega_const_1", "samedata"), sp("1_EFT_GR", "samedata"),
                 sp("5_quint_exponential", "samedata"), sp("4_Des_F_R_LCDM_1", "samedata")],
    # one CAMBparams reused, only an EFT parameter value changed through initialize_parameters
    "reinit_value": [sp("2_PEFT_Omega_const_1"), sp("2_PEFT_Omega_const_1", "reinit", EFTOmega0=0.2),
                     sp("2_PEFT_Omega_const_1", "reinit"), sp("5_quint_exponential", "reinit"),
                     sp("2_PEFT_Omega_const_1", "reinit", EFTOmega0=0.2)],
    # copies of a CAMBparams carrying EFTCAMB (each in its own process: may crash)
    "copy_fortran": [sp("4_Des_F_R_LCDM_1", "copy"), sp("2_PEFT_Omega_const_1", "copy")],
    "deepcopy_peft": [sp("2_PEFT_Omega_const_1", "deepcopy")],
    "deepcopy_fr": [sp("4_Des_F_R_LCDM_1", "deepcopy")],
    "pickle_peft": [sp("2_PEFT_Omega_const_1", "pickle")],
}


def step_key(step: dict) -> str:
    """Identity of a step's physics: steps with the same key must give bit-identical results."""
    if step["kind"] == "ini":
        ov = ",".join(f"{k}={v}" for k, v in sorted(step["overrides"].items()) if k not in FAST)
        return f"ini:{step['model']}:{ov}"
    nonu = "nonu" if step["cosmo"].get("mnu", 0.06) == 0 else "nu"
    eft = ",".join(f"{k}={v}" for k, v in sorted(step.get("eft", {}).items()))
    return f"sp:{step['model']}:{nonu}" + (f":{eft}" if eft else "")


def ini_path(model: str) -> Path:
    p = harness.PARAMS_DIR / f"{model}.ini"
    return p if p.exists() else HERE.parents[1] / "inis" / f"{model}.ini"


def eft_keys(model: str) -> dict[str, str]:
    """EFT keys of a model: all keys of the model ini itself (not of base_params), minus output settings."""
    keys = {}
    for raw in ini_path(model).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" in line and not line.upper().startswith(("DEFAULT", "INCLUDE")):
            k, v = (s.strip() for s in line.split("=", 1))
            keys[k] = v
    for k in ("output_root", "output_file_headers", "EFTCAMB_write_background"):
        keys.pop(k, None)
    keys.setdefault("EFTflag", "0")
    keys["feedback_level"] = 0  # int: set_params compares it with 0
    return keys


def fingerprint(camb, cp, results) -> dict[str, np.ndarray]:
    out = {}
    cls = results.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)
    for key in ("total", "unlensed_scalar", "tensor", "lens_potential"):
        if key in cls:
            out[f"cl_{key}"] = np.asarray(cls[key])
    if cp.WantTransfer:
        k, z, pk = results.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
        out["pk_k"], out["pk"] = np.asarray(k), np.asarray(pk)
        out["sigma8"] = np.atleast_1d(results.get_sigma8())
    out["bg_H"] = np.asarray(results.hubble_parameter(harness.Z_BG))
    out["bg_tau"] = np.asarray(results.conformal_time(harness.Z_BG))
    der = results.get_derived_params()
    out["derived"] = np.array([der[k] for k in sorted(der)])
    if cp.EFTCAMB.EFTflag != 0:
        _, arr = cp.EFTCAMB.get_eft_functions(results, A_GRID)
        for f in arr.dtype.names:
            out[f"eft_{f}"] = np.ascontiguousarray(arr[f])
    return out


def build_setparams(camb, step: dict):
    params = dict(step["cosmo"])
    params.update(eft_keys(step["model"]))
    params.update(step.get("eft", {}))
    dropped = []
    for _ in range(50):
        try:
            return camb.set_params(**params), dropped
        except camb.baseconfig.CAMBUnknownArgumentError as e:
            key = str(e).split(":")[-1].strip()
            if key not in params:
                raise
            dropped.append(key)
            params.pop(key)
    raise RuntimeError("too many unknown keys")


def run_sequence(name: str, steps: list[dict], outdir: Path) -> None:
    """Runs inside the subprocess."""
    sys.path.insert(0, str(harness.tree_path(TREE)))
    import copy
    import pickle

    import camb

    outdir.mkdir(parents=True, exist_ok=True)
    cp = None
    shared_data = None
    log = []
    for i, step in enumerate(steps):
        t0 = time.time()
        rec = {"i": i, "key": step_key(step), "kind": step["kind"], "status": "ok"}
        try:
            kind = step["kind"]
            if kind == "ini":
                wrapper = harness._wrapper_ini(ini_path(step["model"]), outdir, f"{name}_{i}", step["overrides"])
                cp = camb.read_ini(str(wrapper))
                calc = cp
            elif kind == "setparams":
                cp, rec["dropped"] = build_setparams(camb, step)
                calc = cp
            elif kind == "samedata":
                cp, rec["dropped"] = build_setparams(camb, step)
                calc = cp
            elif kind == "samecp":
                cp, rec["dropped"] = build_setparams(camb, step)
                camb.get_results(cp)
                calc = cp
            elif kind == "reinit":
                keys = eft_keys(step["model"])
                keys.update(step.get("eft", {}))
                cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
                calc = cp
            elif kind in ("deepcopy", "copy", "pickle"):
                cp, rec["dropped"] = build_setparams(camb, step)
                if kind == "deepcopy":
                    calc = copy.deepcopy(cp)
                elif kind == "copy":
                    calc = cp.copy()
                else:
                    calc = pickle.loads(pickle.dumps(cp))
                rec["copy_EFTflag"] = int(calc.EFTCAMB.EFTflag)
                rec["copy_num_params"] = int(calc.EFTCAMB.num_params()) if calc.EFTCAMB.EFTflag else 0
            else:
                raise ValueError(kind)
            rec["EFTflag"] = int(calc.EFTCAMB.EFTflag)
            if kind == "samedata":
                if shared_data is None:
                    shared_data = camb.CAMBdata()
                shared_data.calc_power_spectra(calc)
                results = shared_data
            else:
                results = camb.get_results(calc)
            fp = fingerprint(camb, calc, results)
            np.savez(outdir / f"{name}_{i}.npz", **fp)
        except Exception as e:  # noqa: BLE001
            rec["status"] = f"error:{type(e).__name__}:{str(e).strip()[:200]}"
            rec["trace"] = traceback.format_exc()[-800:]
        rec["elapsed"] = round(time.time() - t0, 2)
        log.append(rec)
        print(json.dumps({k: v for k, v in rec.items() if k != "trace"}), flush=True)
        (outdir / f"{name}.log.json").write_text(json.dumps(log, indent=1))


def spawn(name: str, steps: list[dict], outdir: Path) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONHASHSEED="0")
    env.pop("PYTHONPATH", None)
    cfg = json.dumps({"name": name, "steps": steps, "outdir": str(outdir)})
    t0 = time.time()
    with harness.slot(1), open(outdir.parent / f"{outdir.name}_{name}.stdout", "w") as f:
        try:
            rc = subprocess.run([sys.executable, __file__, "_child", cfg], stdout=f, stderr=subprocess.STDOUT,
                                env=env, cwd=outdir, preexec_fn=harness._nice, timeout=540).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return {"name": name, "rc": rc, "elapsed": round(time.time() - t0, 1)}


def fresh_name(key: str) -> str:
    return "fresh_" + "".join(c if c.isalnum() else "_" for c in key)


def compare_fp(pa: Path, pb: Path) -> tuple[bool, float, str]:
    A, B = np.load(pa), np.load(pb)
    exact, worst, worst_k = True, 0.0, ""
    for k in sorted(set(A.files) | set(B.files)):
        if k not in A.files or k not in B.files:
            exact, worst, worst_k = False, float("inf"), f"missing:{k}"
            continue
        x, y = A[k], B[k]
        if x.shape != y.shape:
            exact, worst, worst_k = False, float("inf"), f"shape:{k}"
            continue
        if not np.array_equal(x, y, equal_nan=True):
            exact = False
            scale = np.nanmax(np.abs(y)) if np.any(y) else 1.0
            r = float(np.nanmax(np.abs(x - y)) / scale) if scale else float("inf")
            if r > worst:
                worst, worst_k = r, k
    return exact, worst, worst_k


def main_all(names: list[str]) -> None:
    from concurrent.futures import ThreadPoolExecutor

    WORK.mkdir(parents=True, exist_ok=True)
    seqdir = WORK / TREE / "seq"
    freshdir = WORK / TREE / "fresh"
    seqdir.mkdir(parents=True, exist_ok=True)
    freshdir.mkdir(parents=True, exist_ok=True)
    # fresh single-step runs for every distinct step key
    fresh_steps = {}
    for n in names:
        for st in SEQUENCES[n]:
            key = step_key(st)
            if st["kind"] in ("reinit", "deepcopy", "copy", "pickle", "samecp", "samedata"):
                st = dict(st, kind="setparams")
            fresh_steps.setdefault(key, st)
    jobs = [(fresh_name(k), [st]) for k, st in fresh_steps.items() if not (freshdir / f"{fresh_name(k)}_0.npz").exists()
            and not (freshdir / f"{fresh_name(k)}.log.json").exists()]
    jobs += [(n, SEQUENCES[n]) for n in names]
    if os.environ.get("T2B_REPORT_ONLY"):  # only re-print the comparison table from existing outputs
        jobs = []
    with ThreadPoolExecutor(int(os.environ.get("T2B_JOBS", "4"))) as ex:
        futs = [ex.submit(spawn, n, st, freshdir if n.startswith("fresh_") else seqdir) for n, st in jobs]
        for f in futs:
            r = f.result()
            print(f"  ran {r['name']:50s} rc={r['rc']} {r['elapsed']}s", flush=True)
    summary = {}
    for n in names:
        logp = seqdir / f"{n}.log.json"
        log = json.loads(logp.read_text()) if logp.exists() else []
        first_seen: dict[str, int] = {}
        rows = []
        for i, st in enumerate(SEQUENCES[n]):
            key = step_key(st)
            rec = log[i] if i < len(log) else {"status": "not-run (process died?)"}
            row = {"i": i, "key": key, "kind": st["kind"], "status": rec["status"]}
            for extra in ("copy_EFTflag", "copy_num_params", "dropped"):
                if extra in rec:
                    row[extra] = rec[extra]
            p = seqdir / f"{n}_{i}.npz"
            fresh_p = freshdir / f"{fresh_name(key)}_0.npz"
            fresh_log = freshdir / f"{fresh_name(key)}.log.json"
            fresh_status = json.loads(fresh_log.read_text())[0]["status"] if fresh_log.exists() else "missing"
            row["fresh_status"] = fresh_status
            if p.exists() and fresh_p.exists():
                row["vs_fresh"] = compare_fp(p, fresh_p)
            elif rec["status"].startswith("error") != fresh_status.startswith("error"):
                row["vs_fresh"] = (False, float("inf"), "status differs from fresh run")
            if key in first_seen and p.exists():
                q = seqdir / f"{n}_{first_seen[key]}.npz"
                if q.exists():
                    row["vs_first"] = compare_fp(p, q)
            first_seen.setdefault(key, i)
            rows.append(row)
            vf = row.get("vs_fresh")
            v1 = row.get("vs_first")
            print(f"{n:18s} {i} {st['kind']:9s} {key[:60]:60s} {rec['status'][:40]:40s} "
                  f"fresh={'-' if vf is None else ('EXACT' if vf[0] else f'{vf[1]:.2e}:{vf[2]}')} "
                  f"first={'-' if v1 is None else ('EXACT' if v1[0] else f'{v1[1]:.2e}:{v1[2]}')}")
        summary[n] = rows
    (WORK / TREE / f"summary_{'_'.join(names) if len(names) < 3 else 'all'}.json").write_text(
        json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    if sys.argv[1] == "_child":
        harness._exit_with_parent()
        cfg = json.loads(sys.argv[2])
        run_sequence(cfg["name"], cfg["steps"], Path(cfg["outdir"]))
    elif sys.argv[1] == "all":
        main_all(sys.argv[2].split(",") if len(sys.argv) > 2 else list(SEQUENCES))
    elif sys.argv[1] == "seq":
        main_all([sys.argv[2]])
