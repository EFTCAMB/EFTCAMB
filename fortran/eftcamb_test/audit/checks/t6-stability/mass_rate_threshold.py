"""(iv) EFT_mass_stability_rate and EFTCAMB_stability_threshold semantics."""
import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
from scan_all import sampler_grid
from boundary_scans import bisect, code_stable
ini = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
lin, log1, _ = sampler_grid()
g = np.unique(np.concatenate([lin, log1]))
r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=g, feedback=0,
                cleanup_each=False)
a = r["a"]; m = a > r["pert_turn_on"]
x = np.minimum(r["EFT_mu1"], r["EFT_mu2"]) * (a / r["adotoa"]) ** 2
expected = float(np.sqrt(-x[m].min()))
rate, slo, shi, n = bisect("A0", ini, "EFT_mass_stability_rate", 0.5, 3.0, {}, rtol=1e-5)
out = {"model": ini.stem, "expected_rate_from_sampled_mu": expected, "code_rate_boundary": rate,
       "rel": abs(rate - expected) / expected, "verdict_lo_hi": [slo, shi], "a_min": float(a[m][np.argmin(x[m])])}
# threshold: does EFTCAMB_stability_threshold change any verdict?
unst = {"EFTGamma10": "-1.0"}
s0, e0 = code_stable("A0", ini, unst)
s1, e1 = code_stable("A0", ini, {**unst, "EFTCAMB_stability_threshold": "1e30"})
s2, e2 = code_stable("A0", ini, {**unst, "EFTCAMB_stability_threshold": "-1e30"})
out["threshold"] = {"default": s0, "1e30": s1, "-1e30": s2}
print(json.dumps(out, indent=1))
