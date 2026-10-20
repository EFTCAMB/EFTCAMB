"""Independent verification of finding t3-counts-rsd-gr-branch-missing-de.

Differences from the t3-gr repro:
  * different EFT model: the shipped test model 2_PEFT_wDE_LCDM_1 (Pure EFT, EFTwDE=1 with w0=-1, Omega=0) and the
    shipped 2_PEFT_gamma4_lin_1 (RGR time 1.1 in the A0 reference log, i.e. EFT perturbations never on);
  * different reference: the same tree with EFTflag=0 (1_EFT_GR.ini), so no B0-vs-A0 differences enter at all;
  * the redshift-distortion source is isolated (all other number-count source terms and counts lensing off);
  * the pre-rebase tree B1 is checked the same way (EFTCAMB initialised from the flattened ini, as the harness does);
  * analytic check: from the EFT timestep cache (get_eft_functions) compare the GR-branch coefficient
    4 adotoa^2 + gpres + grho/3 (gpres = gpresm_t, grho = grhom_t, the values derivs passes to output_window_sources
    when EFTflag/=0) with the Friedmann/Raychaudhuri value 4 adotoa^2 - 2 Hdot, around the window redshift.

Usage: python t3counts_indep.py [outdir]
"""

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

PAR = harness.PARAMS_DIR
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "_worker" else (
    harness.AUDIT_ROOT / "work" / "verify" / "t3counts" / "indep")

RSD_ONLY = {"counts_density": "F", "counts_radial": "F", "counts_timedelay": "F", "counts_ISW": "F",
            "counts_velocity": "F", "counts_potential": "F", "DoRedshiftLensing": "F", "counts_redshift": "T"}
NO_RSD = {"counts_redshift": "F"}


def worker(tree, ini, out):
    sys.path.insert(0, tree)
    import camb

    cp = camb.read_ini(ini, no_validate=True)
    info = {"camb_file": camb.__file__}
    keys = harness.parse_ini(ini)
    if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag == 0 and int(keys.get("EFTflag", "0")) != 0:
        cp.EFTCAMB.initialize_parameters(cp, keys, print_header=False)
        info["eft_init_fallback"] = True
    info["EFTflag"] = int(cp.EFTCAMB.EFTflag) if hasattr(cp, "EFTCAMB") else None
    info["max_redshift"] = float(cp.Reion.max_redshift) if hasattr(cp.Reion, "max_redshift") else None
    info["counts_redshift"] = bool(cp.SourceTerms.counts_redshift)
    info["counts_density"] = bool(cp.SourceTerms.counts_density)
    info["nwin"] = len(cp.SourceWindows)
    r = camb.get_results(cp)
    d = r.get_source_cls_dict(raw_cl=True, lmax=800)
    arrs = {k: np.asarray(v) for k, v in d.items()}
    if info["EFTflag"]:
        z = np.array([2.0, 1.0, 0.7, 0.5, 0.3, 0.0])
        a = 1 / (1 + z)
        try:
            fields, arr = cp.EFTCAMB.get_eft_functions(r, a)
            arrs["eft_a"] = a
            arrs["eft_fields"] = np.array(fields)
            arrs["eft_arr"] = np.asarray(arr)
        except Exception as e:  # pre-rebase python bridge may differ
            info["eft_functions_error"] = repr(e)
    np.savez(out, **arrs)
    Path(str(out) + ".json").write_text(json.dumps(info))


def run(name, tree, base, ov):
    OUT.mkdir(parents=True, exist_ok=True)
    ini = OUT / f"{name}.ini"
    ini.write_text("\n".join([f"DEFAULT({base})", f"output_root = {OUT / name}", "EFTCAMB_write_background = F",
                              "feedback_level = 0", "get_tensor_cls = F"] + [f"{k} = {v}" for k, v in ov.items()]))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    with harness.slot(1):
        p = subprocess.run([sys.executable, __file__, "_worker", str(harness.tree_path(tree)), str(ini),
                            str(OUT / f"{name}.npz")], env=env, capture_output=True, text=True, cwd=OUT)
    if p.returncode:
        print(name, "FAILED", p.stderr[-1500:])
    return name, p.returncode


def rel(a, b, k):
    return float(np.max(np.abs(a[k][2:] - b[k][2:])) / np.max(np.abs(b[k][2:])))


def main():
    gr = PAR / "1_EFT_GR.ini"
    wde = PAR / "2_PEFT_wDE_LCDM_1.ini"
    g4 = PAR / "2_PEFT_gamma4_lin_1.ini"
    cases = []
    for tree in ("A0", "B1"):
        cases += [(f"{tree}_GR", tree, gr, {}), (f"{tree}_wDE", tree, wde, {}),
                  (f"{tree}_GR_rsd", tree, gr, RSD_ONLY), (f"{tree}_wDE_rsd", tree, wde, RSD_ONLY),
                  (f"{tree}_GR_norsd", tree, gr, NO_RSD), (f"{tree}_wDE_norsd", tree, wde, NO_RSD)]
    cases += [("A0_g4lin", "A0", g4, {}), ("A0_wDE_rsd_forced", "A0", wde,
                                            dict(RSD_ONLY, EFTCAMB_GR_threshold="-1", EFTCAMB_skip_stability="T"))]
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda c: run(*c), cases))

    def ld(n):
        p = OUT / f"{n}.npz"
        if not p.exists():
            return None
        with np.load(p) as f:
            d = {k: f[k] for k in f.files}
        d["info"] = json.loads(Path(str(p) + ".json").read_text())
        return d

    results = {}
    pairs = [("A0_wDE", "A0_GR"), ("A0_wDE_norsd", "A0_GR_norsd"), ("A0_wDE_rsd", "A0_GR_rsd"),
             ("A0_wDE_rsd_forced", "A0_GR_rsd"), ("A0_g4lin", "A0_GR"),
             ("B1_wDE", "B1_GR"), ("B1_wDE_norsd", "B1_GR_norsd"), ("B1_wDE_rsd", "B1_GR_rsd")]
    for n, ref in pairs:
        a, b = ld(n), ld(ref)
        if a is None or b is None:
            print(f"{n} vs {ref}: missing")
            continue
        row = {k: rel(a, b, k) for k in ("W1xW1", "W1xW2", "W2xW2", "PxW1")}
        results[f"{n}_vs_{ref}"] = row
        print(f"{n:20s} vs {ref:12s}: " + "  ".join(f"{k}={v:.2e}" for k, v in row.items()),
              "| EFTflag", a["info"]["EFTflag"], b["info"]["EFTflag"], "fallback", a["info"].get("eft_init_fallback"))

    # analytic coefficient check from the EFT cache
    d = ld("A0_wDE")
    if d is not None:
        f = list(d["eft_fields"])
        arr = d["eft_arr"]
        ix = {k: f.index(k) if k in f else -1 for k in ("a", "adotoa", "Hdot", "grhom_t", "gpresm_t", "grhov_t", "gpiv_t")}
        print("analytic GR-branch coefficient (cache of A0 2_PEFT_wDE_LCDM_1):")
        print("   z     buggy=4H^2+gpres+grho/3   correct=4H^2-2Hdot   rel.err   missing/(-(2/3)grhov_t)")
        coef = []
        for i, a in enumerate(d["eft_a"]):
            v = {k: (arr[i][k] if arr.dtype.names else arr[i, j]) for k, j in ix.items()}
            buggy = 4 * v["adotoa"] ** 2 + v["gpresm_t"] + v["grhom_t"] / 3
            corr = 4 * v["adotoa"] ** 2 - 2 * v["Hdot"]
            miss = corr - buggy
            ratio = miss / (-(2 / 3) * v["grhov_t"]) if v["grhov_t"] else float("nan")
            coef.append({"z": 1 / a - 1, "buggy": buggy, "correct": corr, "rel": buggy / corr - 1, "ratio": ratio})
            print(f"  {1 / a - 1:4.2f}   {buggy:.6e}            {corr:.6e}      {buggy / corr - 1:+.4f}   {ratio:.6f}")
        results["coef"] = coef
    (OUT / "summary.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        worker(*sys.argv[2:5])
    else:
        main()
