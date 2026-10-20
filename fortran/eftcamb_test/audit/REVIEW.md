# EFTCAMB audit: code bugs for review

Code under audit: branch `eftcamb-audit` (f6f4ba9; code identical to 5bef49a). Tracks 0-6 and 8-10 are complete; the
physics workflow (track 7) has not run. IDs refer to `FINDINGS.md` (regenerate with `checks/compile_findings.py`).

All 108 verifier runs confirmed their finding; none was refuted. The 143 ledger entries merge to the items below.
Each S1/S2 item has a reproducer in `repro/`. "Tested" means a track agent applied the fix in a tree copy and
measured before/after; those copies were lost in the scratchpad wipe, so every fix will be re-applied and re-tested.

Severity after verification is shown as S1/S2/S3. For each item: what goes wrong, where, the proposed fix.

## A. Safe-local fixes (crashes, Python interface, uninitialised data)

Expected effect on shipped-ini spectra: none (<= 1e-5); they turn crashes or silent garbage into correct results or
clean errors.

| # | Sev | Problem | Where | Proposed fix | IDs |
|---|---|---|---|---|---|
| A1 | S1 | Reusing `results.Params` (incl. `get_BAO`, which v2 calls on `Params.copy()`) aborts Python: allocate of already-allocated parameter-cache arrays | `fortran/results.f90:580` | Deallocate (or skip) before allocating | F-001, F-004 |
| A2 | S1 | `pickle`/`deepcopy`/`repr`+`eval` of CAMBparams drop the Fortran EFT model; next run segfaults (B1 too) | `camb/model.py`, `camb/baseconfig.py:717` | `__getstate__` stores `read_parameters()`, `__setstate__` re-initialises EFTCAMB | F-002, F-005 |
| A3 | S1 | Python getters segfault when no model/cache exists: `num_params`/`param_names` with EFTflag=0, after a failed `initialize_parameters`, `get_eft_functions` on GR results, `get_scale_evolution` after a transfer-only run | `fortran/eftcamb/11_EFTCAMB_python.f90:200,325,382` | Guard every getter with `allocated(...)`; reset EFTflag=0 on failed init; force WantScalars inside GetEvolution | F-007, F-021, F-024 |
| A4 | S1 | 3+1 wDE designer always fails ("background solver failed"): intent(out) `success` never set | `07f_designer_models/007p8_3p1_wDE_designer.f90:888` | `success = .true.` at entry (one line) | F-008 |
| A5 | S2 | `get_evolution`/`get_scale_evolution` sort the requested times silently; rows misaligned for decreasing a | `camb/eftcamb.py:316` | argsort in, inverse argsort out | F-022, F-081 |
| A6 | S2 | `set_params(cp=existing, ...)` silently resets the EFT model to GR; EFT-only update raises | `camb/camb.py:216-236` | Merge existing `read_parameters()` with the new keys; skip when no EFT key given | F-045, F-079 |
| A7 | S2 | `write_ini` writes no EFTCAMB keys: validation fails for every EFT model, `validate=False` writes a GR ini | `camb/_ini.py:106` | Write `read_parameters(raw=True)` | F-080, F-089 |
| A8 | S2 | `read_parameters` returns bytes for strings and Fortran `1.d-4` numbers; feeding the dict back garbles them | `camb/eftcamb.py:184` | Decode, handle D exponents, add `raw=True` | F-078 |
| A9 | S2 | Fortran `stop`/`MpiStop`/`error stop` kill the Python interpreter: 131 sites (18 reachable per k), unparsable ini values, documented `EFTwDE` 6/7/11 in designer models | many; e.g. `008p3_Kmouflage.f90:438`, `11_EFTCAMB_python.f90:36`, `007p3_Designer_fR.f90:471` | Route through `GlobalError`/`eft_error`; Python-side ini reader that sets a flag instead of stopping. **Large** (touches ~26 files; see D10) | F-036, F-077, F-047 |
| A10 | S2 | Equispaced linear interpolation reads index N+1 when x is within rounding of x_final (hit by `3_RPH_KT_power_1` under bounds checking) | `02_equispaced_interpolation_linear_1D.f90:234` (7 places) | Clamp the index | F-025 |
| A11 | S2 | Command-line camb segfaults in debug builds; `SetParams` tests the previous parameters' EFTflag | `fortran/results.f90:365` | Test `P%EFTCAMB` instead of `this%CP%EFTCAMB` | F-026, F-090 |
| A12 | S2 | Timestep-cache init leaves `EFTcdotdot`, `EFTGamma1PP` (Horndeski `V_scf`) unset: values returned to Python and NaN-checked by the stability code are garbage | `03_abstract_EFTCAMB_cache.f90:402,443` | Zero them in `TimestepCacheInit` | F-031, F-065 |
| A13 | S2 | Stability routine never resets `gpinu_tot` in its 2nd/3rd neutrino loops: tripled (wrong Hdotdot inputs to stability decisions only) | `09_EFTCAMB_stability.f90:493,522` | Reset `gpinu_tot`, not `gpinudot_tot` | F-066 |
| A14 | S2 | `EFT_IC_type` 2/3 read Gamma1-3 and background from a stale or uninitialised cache; type 3 returns all-NaN spectra with no error when Omega=Gamma3=0 | `fortran/equations.f90:2243`, `09_EFTCAMB_IC.f90:325` | Fill the cache before use; guard the zero denominator; NaN check after IC | F-056, F-057 |
| A15 | S2 | `EFTCAMB_turn_on_time` = 0, negative, NaN or >1 silently runs GR perturbations | `09_EFTCAMB_RGR.f90:143`, `09_EFTCAMB_main.f90:240` | Validate in (0,1], raise | F-082 |
| A16 | S2 | Early-turn-on models (5_ADE_1, a=1e-10) refused with a hard error when `l_max_scalar` is lowered | `fortran/equations.f90:2188` | Start thermal history at min(tauminn, turn-on) | F-044 |
| A17 | S2 | theta->H0 solve aborts when an EFT model is unstable at a bracket end | `camb/model.py:786` | Shrink the failing bracket end before brentq | F-071 |
| A18 | S2 | `import camb` fails for ifort/flang libraries (wrong bind(c) symbol name) | `camb/eftcamb.py:23` | Use the bind(c) name unconditionally | F-084 |
| A19 | S3 | Horndeski `integer :: num_eq = 5` is an implicit SAVE: 6 leaks into later runs | `008p0_Horndeski.f90:2381` | Separate declaration and assignment | F-102, F-105 |
| A20 | S3 | Uninitialised locals in EFT blocks: `grhog_t`/`grhor_t` in `derivst`, `w_dark_energy_t` in `derivs` | `fortran/equations.f90:3688,4263` | Assign before use | F-092, F-093, F-098, F-129 |
| A21 | S3 | String getters called without gfortran's hidden length argument (works only because Python allocates exactly 200 bytes); immutable bytes used as output buffers | `11_EFTCAMB_python.f90:248`, `camb/eftcamb.py:211,267` | bind(c) wrappers with explicit length, `create_string_buffer` | F-085, F-087, F-140 |
| A22 | S3 | With EFTflag=0, `set_params(EFTflag=0, <EFT key>)` raises CAMBUnknownArgumentError (user report) | `camb/camb.py:224` | Mark known EFT keys as used | F-142 |
| A23 | S3 | Smaller items: RPH getter sizes a block with the wrong count (F-104); QS outputs 0/0 NaN (F-103); OpenMP-private EV flags garbage in debug path (F-091, F-101); 2D feedback writes absent optional (F-118); upstream `size()` of unallocated array (F-086) | various | One-line guards | F-086, F-091, F-101, F-103, F-104, F-118 |

## B. Physics-changing fixes with a clear correct answer

Typos and missing terms where the right expression is unambiguous; each changes outputs of some models, so each needs
your sign-off. Size of the change is from the reproducers.

| # | Sev | Problem | Where | Fix | Effect | IDs |
|---|---|---|---|---|---|---|
| B1 | S1 | `EFTpiDfunction = D1 + k^2 D1` instead of `D1 + k^2 D2`: wrong initial pi' | `09_EFTCAMB_IC.f90:589` | one character | K-mimic P(k) at k=1-20/Mpc changes by 0.3-1.3 (e.g. +1.36 -> +0.06 at k=10) | F-015, F-113 |
| B2 | S1 | `EFT_kinetic` multiplies the `6 H gamma2 Omega'` cross term by h0^2 instead of h0: scalar ghost check wrong for any model with gamma2 and Omega' (RPH: up to 6x off; a ghost model accepted) | `06_abstract_EFTCAMB_model.f90:690` | move a parenthesis | Stability verdicts; all 73 shipped verdicts unchanged | F-017 |
| B3 | S1 | `derivs` uses pi'' in `EFTeomM`/`EFTeomNdot` before computing it: ISW source, sigma'', Psi' miss the pi'' term | `fortran/equations.f90:2989` | recompute the two factors after the pi equation | Horava TT up to 4.6e-3 (ell 4), SolSyst_3 -2.6% at ell 2-3 | F-018 |
| B4 | S2 | Third stability sampling pass never updates its sample point: re-tests a=astart 1000 times; a ghost just below a=1 is accepted | `09_EFTCAMB_stability.f90:127` | add the missing `y = ...` line | Stability verdicts only | F-035, F-061 |
| B5 | S2 | QS Sigma output is a verbatim copy of the QS mu formula | `06_abstract_EFTCAMB_model.f90:972` | implement Sigma | `get_evolution` output only | F-034 |
| B6 | S2 | Full-mapping `Hdotdotdot` wrong (typos, wrong prefactor, uses EFTLambdadotdot before it is computed; 223-822x off for K-mouflage); designer f(R) leaves it 0; other higher derivatives unset | `06p1_abstract_EFTCAMB_full_map.f90:143-168` | exact derivative of the Hdotdot expression | Positivity bounds and Python output; not Cls | F-023, F-027, F-054 |
| B7 | S1 | Hořava `compute_H_derivs`: Hdot/Hdotdot inconsistent with the model's own adotoa (4-62% off; second Friedmann residual up to 25%) | `008p1_Horava.f90:453` | rewrite from the model's H(a) | Hořava spectra | F-013 |
| B8 | S2 | Scaling Cubic third derivatives: hpprime sign and y2 coefficient, missing x1pp term | `008p6_Scaling_Cubic.f90:931` | per diff in finding | 5_scg higher derivatives | F-053 |
| B9 | S2 | EFT rewrite of `polterddot` drops the curvature factor `Kf(2)`: EFTflag=0 non-flat CMB differs from CAMB (TT 6.4e-3 at ell 2, omk=-0.05) | `fortran/equations.f90:3858` | restore `Kf(2)` | GR non-flat matches CAMB to 1e-13 | F-037, F-068 |
| B10 | S2 | Number-counts RSD source misses dark energy for EFT models in the GR branch (6.9e-3 on counts Cls) | `fortran/equations.f90:1698` | use `-2 Hdot` identity | counts Cls -> 1e-7 of CAMB | F-038 |
| B11 | S2 | Stability check evaluates dtauda with zero matter density (`grhoa2` never set): near-GR models rejected with "dtauda is NaN"; `grhoa2` also wrong in `derivs` and the Python bridge | `09_EFTCAMB_stability.f90:191`, `equations.f90:2792`, `11_EFTCAMB_python.f90:418` | set `grhoa2 = grhom_t a^2` | Verdicts, Python output | F-043 |
| B12 | S2 | `compute_additional_derivs` differentiates helpers at the host a/cache, so mass stability sees Gamma1PP=Gamma2PP=0 etc.; Minkowski-limit/positivity checks run before it | `06p2_abstract_EFTCAMB_designer.f90:141`, `09_EFTCAMB_stability.f90:570` | pass the temporary a/cache; call it first | Verdicts of optional checks | F-032, F-062 |
| B13 | S1/S2 | Parametrization maths: cubic spline wrong at end nodes and derivatives inconsistent (F-011); splines return `null_value` instead of 0 outside nodes (F-048); quintic spline `dx**6` instead of `dx**5` (F-049); spline integrals NaN/wrong so designer f(R) cannot use spline w (F-046); quintic order 3 ERROR STOP (F-050); steplog d4 wrong/NaN, integral(1)!=1 (F-012); Padé d4 drops a factor (F-052); latent d4 errors in unused types (F-119) | `04f_parametrizations_1D/*`, `02p12_splines5.f90:391` | per finding (all tested to machine precision) | Only runs using these parametrizations | F-011, F-012, F-046, F-048, F-049, F-050, F-052, F-119 |
| B14 | S2 | Ini keys truncated to 20 characters: documented long keys (e.g. `Designer3p1sigma1_ODE0`, `RPHbraiding_ODE_delta`) silently read as 0 | `04_abstract_parametrizations_1D.f90:343` | raise `EFT_names_max_length` | Any ini using long keys reads different values | F-051 |
| B15 | S1 | Massive neutrinos missing from model backgrounds: Extended Galileon (for a>=0.2; also no DE in CAMB H(z) at z>4, NaN at S=0), Scaling Cubic (H(a) off 3.5-7% vs CAMB) | `008p7_Extended_Galileon.f90:435,455`, `008p6_Scaling_Cubic.f90:1123` | compute nu density locally / add to the ODE, or refuse massive nu | These two models | F-006, F-010, F-014 |
| B16 | S2 | Full-mapping `derivs` evaluates background EFT functions before adotoa/Hdot/Hdotdot are set (Beyond Horndeski, Extended Galileon read stale values) | `fortran/equations.f90:2816` | compute from the model's own background | Those two models | F-033 |
| B17 | S1 | thetastar/cosmomc_theta ini keys (new in v2) solve H0 before EFTCAMB is initialised: wrong H0 for every EFT model using theta in an ini | `fortran/camb.f90:711` | initialise EFTCAMB first | Ini runs with theta | F-020 |
| B18 | S2 | `get_background_densities`, `get_Omega('de')`, `get_dark_energy_rho_w` ignore the EFT background (H rebuilt from them off by up to 5%) | `fortran/results.f90:1157` | use EFT background | Those outputs | F-070 |
| B19 | S2 | EFTCAMB's `get_matter_power_interpolator(extrap_kmax)` change extrapolates ~10x worse than upstream (+3.4 vs +0.35 at 10 kmax) | `camb/results.py:1188` | revert to upstream slope (+ non-log branch fix) | Extrapolated P(k) | F-069 |
| B20 | S2 | `EFTCAMB_use_background=T`: cosmomc_theta uses the exact sound horizon instead of CAMB's cosmomc approximation (shift -1e-3, ~4 sigma Planck) | `fortran/results.f90:1135` | drop that branch | Horndeski theta | F-072 |
| B21 | S2 | K-mouflage/K-mimic, ADE, quintessence background tables too coarse: ~1e-4 Cls, 3-7e-5 H(z) error floor | `008p3_Kmouflage.f90` etc. | ~10x points or cubic interpolation | 1e-4 level | F-041 |
| B22 | S2 | EFT_IC_type 2/3 apply a radiation-era super-horizon series at any turn-on time: garbage spectra (TT 7e5, sigma8 0.55->2.76) for late turn-on | `09_EFTCAMB_IC.f90:307` | apply only when k tau << 1 in radiation era, else type 1 or error | IC 2/3 runs | F-055 |
| B23 | S3 | Kind-less single-precision literal in Scaling Cubic changes results above 1e-8 | `008p6_Scaling_Cubic.f90:247` | add `_dl` | 5_scg at 1e-8 | F-106 |

## Decisions (user review, 2026-10-20)

**Status: all items approved.** Groups A, B and D are accepted as proposed in the tables below; group C is
resolved as listed here. Nothing has been implemented yet.

| Item | Decision | Implementation plan |
|---|---|---|
| A5 | EFTCAMB bug (`camb/eftcamb.py` is EFTCAMB-only; upstream `get_time_evolution` already argsorts and restores order) | Follow upstream: argsort in, inverse argsort out |
| A9 / C12 | Hard inconsistencies must stop the calculation with an error, but must not kill Python | Convert all 131 sites (scope: full). Each `stop`/`MpiStop`/`error stop` becomes `call GlobalError(msg, error_...)` + `return`; callers check `global_error_flag` and return; per-k OpenMP loops already skip work once the flag is set. Python raises `CAMBError(msg)` (existing CAMB mechanism); the command-line program still exits with the message. Ini parsing: a `TIniFile` subclass whose `FailStop` sets the error instead of stopping, which requires `type(TIniFile)` -> `class(TIniFile)` in the EFT dummy arguments (~26 files) |
| B7 | Confirmed against arXiv:1508.01787: code H(a) = paper eq. (Hubble) exactly (F = G_cosmo/G_N = (2+2xi-eta)/(3 lambda+2) with code lambda = lambda_p - 1, xi = xi_p - 1; solar-system eta = 2 xi equals paper eta = 2 xi_p - 2), but the paper's second Friedmann eq. (friedmann2), equivalently d/dtau of eq. (Hubble), gives Hdot = -adotoa^2/2 - F (gpresm_t + gpiv_t)/2 and Hdotdot = -adotoa Hdot - F (d gpres/dtau)/2. The code instead multiplies the GR expressions (already containing the F-scaled adotoa) by sqrt(F): wrong unless F = 1. Unchanged since pre-rebase (B1) | Replace with the derivative of eq. (Hubble); add Hdotdotdot |
| C1 | No curvature | Refuse omk != 0 when EFTflag != 0 with a clear error (results.f90 EFT block) |
| C2 | Yes | Regularise type-1 IC: pi = pi' = 0 when the denominator is small relative to its terms or has no stable attractor; pi' = 0 when coefficients vary fast; k-independent IC time; warn |
| C3 | Yes | Normalise c, Lambda by 3 adotoa^2 (dimensionless RGR test); observable-based criterion for scale-dependent models (designer f(R)) |
| C4 | Recommended option: (a) + (b) + (c) | (a) fix the derivative inputs now (B12); re-derive mu_1, mu_2 in the physics workflow against the reference equations, with a canonical-quintessence unit test (\|mu\|/H^2 = O(1)) and a form that does not divide by vanishing kinetic factors. (b) Until then, `EFT_mass_stability = T` prints a warning that the check is experimental and is ignored. (c) Remove `EFT_mass_stability = T` from the 14 shipped test inis (2_PEFT_gamma1_const_1/2, 2_PEFT_gamma4_lin_1, 2_PEFT_gamma4_power_1/2, 2_PEFT_Omega_const_1/2, 2_PEFT_Omega_lin_1/2, 3_RPH_KT_OmegaDE_1, 3_RPH_KT_power_1, 3_RPH_MB_const_1, 3_RPH_MK_const_1, 3_RPH_MK_power_1) |
| C5 | Disable | `EFTCAMB_evolve_metric_h = T` raises an "unsupported" error |
| C6 | Pass effective w0/wa | Complete `EFTCAMB_effective_w0wa_values` for every model (designer: from the wDE parametrization; full mapping: CPL fit of the background w(a)), feed Halofit/HMcode. Also apply to the v2 automatic non-linear lensing (F-074): keep it, now with the model's w0/wa. 21cm: refuse in the common init (F-075) |
| C7 | Adopt | After `initialize_background`, set CP%H0 from the EFT background and redo the H0-dependent setup |
| C8 | Threshold | New parameter (reuse `EFTCAMB_stability_threshold`, currently read but unused): every stability condition passes if value >= -threshold * (scale of its terms), i.e. a relative numerical floor; default from the measured cancellation level (RPH early-time noise ~1e-10 relative, so ~1e-8) |
| C9 | OK if runtime does not grow much | Implement step limit and turn-on time samples; measure runtime on all 73 models, accept if the median grows < ~10-20% |
| C10 | Intended | Keep; document in the EFTCAMB README |
| C11 | Raise | `set_classes(dark_energy_model!=EFTCAMB)` / `set_dark_energy` with an active EFT model raise CAMBValueError |

## Implementation plan

One commit per item (or per tightly coupled group) on `eftcamb-audit`, never pushed. Order, so that later physics
changes are measured on a crash-free, correctly reporting code base:

1. **Safe-local, group A** (A1-A8, A10-A23), then **A9/C12** (stop -> GlobalError, largest change, its own commits per
   module). Gate: every shipped ini bit-identical to the pre-fix build (or <= 1e-5 where an uninitialised value was
   read), each A reproducer now passes.
2. **Refusals and switches:** C1 (curvature), C5 (metric_h), C11 (set_dark_energy), 21cm refusal, C4 (b)+(c).
   Gate: reproducers raise CAMBError with a clear message; other inis unchanged.
3. **Clear physics fixes, group B**, one at a time (B1-B23, B7 per arXiv:1508.01787). Gate per commit: the
   reproducer gives the expected value; full-model regression lists exactly which models moved and by how much.
4. **Design items:** C2, C3, C6, C7, C8, C9 (C9 with runtime measured on all 73 models; reject if median grows >~20%).
5. **Final validation:** `python -m unittest camb.tests.camb_test`, `pre-commit run --all-files`, full regression of all
   73 test models (+ the 6 audit inis) against the pre-fix build with a per-model table of changes, and the fast
   checks kept as `camb/tests/eftcamb_test.py` (< 2 min).

Then the physics workflow (track 7) runs on the fixed code, including the C4 (a) re-derivation of the mass check.

Resource rules as before: one build at a time, CAMB runs through the harness CPU slots, at most a few agents.

## C. Needs a design decision first (original list; see decisions above)

| # | Sev | Question | Findings | Options suggested by the audit |
|---|---|---|---|---|
| C1 | S1/S2 | **Curvature.** Any EFT model with omk != 0 gives O(1)-wrong spectra (EFT adotoa and perturbation background are flat), silently accepted | F-009, F-076 | (a) refuse omk != 0 when EFTflag != 0 (one block in `results.f90`); (b) implement curvature in all `compute_adotoa`/H derivatives |
| C2 | S1/S2 | **pi initial conditions (type 1).** pi = -E/(C + k^2 D1 + k^4 D2) is singular where the denominator crosses zero (single-k P(k) spikes, e.g. K-mimic factor ~3 at k=0.22); pi' from finite differences of fast-varying coefficients puts pi off its attractor (K-mimic spectra depend on turn-on time by 1e8) | F-003, F-016, F-030, F-040 | regularise (pi=pi'=0 when the denominator is small or no stable attractor), pi'=0 when coefficients vary fast, impose ICs at a k-independent time; tested: glitches vanish, K-mimic converges to GR level |
| C3 | S2 | **Return-to-GR criterion.** Compares dimensionful c, Lambda (Mpc^-2) with a dimensionless threshold; for designer f(R) near GR the threshold, not physics, delays turn-on and moves P(k>=1/Mpc) by 2-11% | F-039, F-058, F-059 | normalise by 3H^2; make it observable-based (e.g. Omega (k_max/aH)^2); or exempt scale-dependent models |
| C4 | S2 | **Mass-stability check.** Rejects canonical quintessence (|mu|/H^2 ~ 1e9) and 51/73 shipped models when forced on; constant Gamma1 -> 0 rejected; positivity allowlist in `camb.py` does not match model names | F-042, F-063, F-124 | re-derive the mass eigenvalues (physics workflow), mark experimental and disable meanwhile, fix allowlist |
| C5 | S2 | **`EFTCAMB_evolve_metric_h`.** Starts from the GR h' and violates the 00 constraint; inconsistent with the constraint formulation (TT 1%, P(k) x13 at k=100) | F-060, F-067 | disable/error until validated, or fix the IC and re-derive the ii-equation path |
| C6 | S2 | **Non-linear and other v2 features with EFT on.** Halofit/HMcode always get w0=-1, wa=0; v2 `set_for_lmax` silently enables GR-calibrated non-linear lensing for MG models; 21cm refused only in the command-line driver | F-073, F-074, F-075 | refuse vs. warn vs. pass effective w0/wa |
| C7 | S1 | **Horndeski H0.** With `Horndeski_shooting=F` the EFT background's own H0 (85.5) never reaches CP%H0 (70.9): sigma8 -12%, h-units wrong | F-019 | adopt the EFT H0 and redo H0-dependent setup, or raise an error when they disagree |
| C8 | S2 | **Ghost-check numerics.** Catastrophic cancellation in `EFT_kinetic` at a < 7e-6 makes a healthy RPH model a ghost at the default stability time; `EFTCAMB_stability_threshold` is read but unused; fixed 3x1000 sampler misses narrow instabilities | F-064, F-125, F-126, F-127, F-128 | cancellation-free alpha-basis form plus a relative tolerance (use the threshold for it); adaptive sampling |
| C9 | S2 | **Integration accuracy of pi.** RK45 steps beyond the stability limit for oscillating pi; EFT turn-on jump (potentials 25-33%) integrated to first order by the CMB time quadrature (low-ell TT 4.6e-2 depends on windows); `get_evolution` pi outputs at k~1/Mpc dominated by integration error | F-028, F-029, F-083 | step limit while EFT active; add time samples straddling the turn-on; tighter tolerance in `get_evolution` |
| C10 | S3 | **Intentional-looking EFTCAMB changes to CAMB.** Reionization `max_redshift` 50 -> 1000 for all runs (only cause of EFTflag=0 vs CAMB difference); EFTflag!=0 forces IntTolBoost>=2; lensing sanity errors turned into silent zero returns; transfer-time mismatch error only when warnings print | F-097, F-111, F-112, F-133, F-134, F-135 | keep and document, or restore upstream |
| C11 | S3 | **`set_dark_energy` / `set_classes(dark_energy_model)` ignored on params with an active EFT model** | F-072 (S2 version), F-132 | raise vs. switch EFT off |
| C12 | S2 | **Scope of the stop -> error refactor (A9).** 131 sites; needs a Python-side ini reader class and `type(TIniFile)` -> `class(TIniFile)` in 26 files | F-036, F-077 | full refactor vs. only the per-k and Python-reachable sites |

## D. Remaining S3 hygiene (not individually listed above)

Grouped hygiene with all locations in `FINDINGS.md`: kind-less literals (F-106), undefined function results and
uninitialised locals (F-107), exact real comparisons in physics branches (F-108), latent thread-safety (F-109),
exact GR points rejected on roundoff (F-110), `want_LCDM` only changes printing (F-114), Fourier/spline5 parameter
edge cases (F-116, F-117), dead and unit-inconsistent EFT code in `equations.f90` (F-130, F-138), massive-nu pressure
perturbation approximation (F-121, F-131), background sampling not doubled for EFT (F-136), Makefile flags (F-137),
H0 search range (F-139), Python input hygiene (F-141), DebugEFTCAMB dump paths (F-096), compiler diagnostics (F-095),
RGR feedback message (F-122), test inis that silently run LCDM (F-123), parameter-cache h0 docs (F-088).
Upstream CAMB only: F-086, F-099, F-100, F-115, F-120.
