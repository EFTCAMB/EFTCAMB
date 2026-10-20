"""Independent verification of nonlinear-halofit-w0wa-lcdm-for-eft (two-sided attribution with an override build).

Uses the instrumented copy built by nlw0wa_build_instr.py: halofit.f90 prints the (w, wa) passed to Halofit / HMcode
on stderr ("AUDIT ...") and EFT_AUDIT_HF_W / EFT_AUDIT_HF_WA override them.

Cases per non-linear model:
  E0  pure-EFT designer DE (no MG), unmodified           -> what EFTCAMB does now
  Ew  same, Halofit/HMcode forced to the true (w0, wa)    -> should equal the fluid
  F   CAMB fluid with the same (w0, wa), unmodified       -> reference
  F1  CAMB fluid, Halofit/HMcode forced to (-1, 0)        -> should equal E0
If |Ew/F - 1| and |E0/F1 - 1| are ~1e-5 while |E0/F - 1| is per cent, the (w, wa) passed to the non-linear model is
the whole cause.  Also: CPL background (w0=-0.9, wa=+0.05), lensed CMB impact, and which (w, wa) other EFT models feed.
"""

import json
import os
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
SCR = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
           "scratchpad/audit")
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else SCR / "work/verify/nlw0wa/instr"
PART = sys.argv[2] if len(sys.argv) > 2 else "all"
OUT = SCR / "work/verify/nlw0wa"
PARAMS = REPO / "fortran/eftcamb_test/parameters"
AINIS = REPO / "fortran/eftcamb_test/audit/inis"
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

sys.path.insert(0, str(TREE))
import camb  # noqa: E402
import numpy as np  # noqa: E402
from camb import model  # noqa: E402

assert Path(camb.__file__).resolve().is_relative_to(TREE.resolve()), camb.__file__
tmp = Path(tempfile.mkdtemp(dir=OUT))
COMMON = "get_tensor_cls = F\nnum_redshiftwindows = 0\nfeedback_level = 0\nEFTCAMB_write_background = F\n"


def capture_stderr(fn):
    """Run fn() with fd 2 redirected to a file; return (result, AUDIT lines)."""
    log = tmp / "stderr.txt"
    sys.stderr.flush()
    saved = os.dup(2)
    with open(log, "w") as f:
        os.dup2(f.fileno(), 2)
        try:
            res = fn()
        finally:
            os.dup2(saved, 2)
            os.close(saved)
    lines = [ln.strip() for ln in log.read_text().splitlines() if "AUDIT" in ln]
    return res, sorted(set(lines))


def run(ini_path, extra, version, override=None, cls=False):
    ini = tmp / "case.ini"
    ini.write_text(f"DEFAULT({ini_path})\n{COMMON}{'get_scalar_cls = T' if cls else 'get_scalar_cls = F'}\n{extra}")
    for key in ("EFT_AUDIT_HF_W", "EFT_AUDIT_HF_WA"):
        os.environ.pop(key, None)
    if override is not None:
        os.environ["EFT_AUDIT_HF_W"], os.environ["EFT_AUDIT_HF_WA"] = (repr(float(v)) for v in override)

    def go():
        cp = camb.read_ini(str(ini))
        cp.set_matter_power(redshifts=[2.0, 1.0, 0.0], kmax=10.0, nonlinear=True)
        cp.NonLinear = model.NonLinear_both if cls else model.NonLinear_pk
        cp.NonLinearModel.set_params(halofit_version=version)
        if cls:
            cp.set_for_lmax(3000, lens_potential_accuracy=1)
        r = camb.get_results(cp)
        k, z, plin = r.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
        _, _, pnl = r.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False)
        out = {"k": k, "z": z, "plin": plin, "ratio": pnl / plin, "eftflag": cp.EFTCAMB.EFTflag}
        if cls:
            out["lensed"] = r.get_lensed_scalar_cls(CMB_unit="muK", raw_cl=True)
        if cp.EFTCAMB.EFTflag != 0:
            out["w0wa_py"] = cp.EFTCAMB.get_effective_w0wa()
        return out

    with harness.slot(1):
        res, audit = capture_stderr(go)
    for key in ("EFT_AUDIT_HF_W", "EFT_AUDIT_HF_WA"):
        os.environ.pop(key, None)
    res["audit"] = audit
    return res


def nl_diff(a, b, kmin=0.05, kmax=5.0):
    sel = (a["k"] > kmin) & (a["k"] < kmax)
    out = {}
    for i, zz in enumerate(a["z"]):
        d = a["ratio"][i, sel] / b["ratio"][i, sel] - 1
        j = int(np.argmax(np.abs(d)))
        out[f"z={zz:.0f}"] = (float(d[j]), float(a["k"][sel][j]))
    lin = float(np.max(np.abs(a["plin"][:, sel] / b["plin"][:, sel] - 1)))
    return out, lin


results = {"tree": str(TREE)}


def pairs_test(tag, eft_ini, eft_extra, fluid_extra, w0, wa, versions):
    for version in versions:
        e0 = run(eft_ini, eft_extra, version)
        ew = run(eft_ini, eft_extra, version, override=(w0, wa))
        f = run(PARAMS / "1_EFT_GR.ini", fluid_extra, version)
        f1 = run(PARAMS / "1_EFT_GR.ini", fluid_extra, version, override=(-1.0, 0.0))
        rec = {}
        for name, (a, b) in {"E0_vs_F (bug)": (e0, f), "Ew_vs_F (EFT forced true w)": (ew, f),
                             "E0_vs_F1 (fluid forced w=-1)": (e0, f1)}.items():
            d, lin = nl_diff(a, b)
            rec[name] = {"nl_ratio_maxdiff": d, "linear_maxdiff": lin}
        rec["audit_E0"], rec["audit_F"] = e0["audit"], f["audit"]
        rec["w0wa_py_E0"] = e0.get("w0wa_py")
        results[f"{tag}:{version}"] = rec
        print(f"[{tag} {version}] audit E0={e0['audit']} F={f['audit']} py w0wa(E0)={e0.get('w0wa_py')}")
        for name in ("E0_vs_F (bug)", "Ew_vs_F (EFT forced true w)", "E0_vs_F1 (fluid forced w=-1)"):
            print(f"    {name:32s} lin={rec[name]['linear_maxdiff']:.2e} nl: "
                  + "  ".join(f"{z}:{v[0]:+.3e}@k={v[1]:.2f}" for z, v in rec[name]["nl_ratio_maxdiff"].items()))
        sys.stdout.flush()


if PART in ("all", "wcdm"):
    pairs_test("wCDM-0.8", PARAMS / "2_PEFT_wDE_wCDM_1.ini", "", "w = -0.8\ncs2_lam = 1\n", -0.8, 0.0,
               ["takahashi", "mead2020", "mead"])

if PART in ("all", "cpl"):
    cpl_eft = "EFTwDE = 2\nEFTw0 = -0.9\nEFTwa = 0.05\n"
    pairs_test("CPL(-0.9,0.05)", PARAMS / "2_PEFT_wDE_wCDM_1.ini", cpl_eft,
               "w = -0.9\nwa = 0.05\ncs2_lam = 1\n", -0.9, 0.05, ["takahashi", "mead2020"])

if PART in ("all", "cmb"):
    e0 = run(PARAMS / "2_PEFT_wDE_wCDM_1.ini", "", "mead2020", cls=True)
    f = run(PARAMS / "1_EFT_GR.ini", "w = -0.8\ncs2_lam = 1\n", "mead2020", cls=True)
    ell = np.arange(e0["lensed"].shape[0])
    rec = {}
    for lmax in (2000, 3000):
        sel = (ell >= 2) & (ell <= lmax)
        a, b = e0["lensed"][sel], f["lensed"][sel]
        rec[lmax] = {s: float(np.max(np.abs(a[:, i] / b[:, i] - 1))) for i, s in enumerate(["TT", "EE", "BB"])}
        rec[lmax]["TE/sqrt(TTEE)"] = float(np.max(np.abs(a[:, 3] - b[:, 3]) / np.sqrt(b[:, 0] * b[:, 1])))
    results["lensedCMB_mead2020_wCDM-0.8_E0_vs_F"] = rec
    print("[lensed CMB, mead2020, NonLinear_both] E0 vs fluid max rel diff:", json.dumps(rec))

if PART in ("all", "scope"):
    scope = {}
    cases = [
        ("4_Des_F_R_wCDM_1", PARAMS / "4_Des_F_R_wCDM_1.ini", ""),
        ("4_Des_mc5e_CPL_1", PARAMS / "4_Des_mc5e_CPL_1.ini", ""),
        ("3_RPH_MK_OmegaDE_1", PARAMS / "3_RPH_MK_OmegaDE_1.ini", ""),
        ("5_quint_exponential", PARAMS / "5_quint_exponential.ini", ""),
        ("5_quint_exponential+effw0wa=T", PARAMS / "5_quint_exponential.ini", "EFTCAMB_effective_w0wa = T\n"),
        ("5_Kmouflage_1", PARAMS / "5_Kmouflage_1.ini", ""),
        ("5_hdsk_TG_1", PARAMS / "5_hdsk_TG_1.ini", ""),
        ("5_hdsk_TG_1+effw0wa=T", PARAMS / "5_hdsk_TG_1.ini", "EFTCAMB_effective_w0wa = T\n"),
    ]
    for name, ini, extra in cases:
        try:
            r = run(ini, extra, "takahashi")
            scope[name] = {"audit": r["audit"], "w0wa_py": r.get("w0wa_py")}
        except Exception as exc:  # noqa: BLE001
            scope[name] = {"error": repr(exc)}
        print(f"[scope] {name}: {scope[name]}")
        sys.stdout.flush()
    results["scope"] = scope

(OUT / f"nlw0wa_indep_{PART}.json").write_text(json.dumps(results, indent=1, default=str))
