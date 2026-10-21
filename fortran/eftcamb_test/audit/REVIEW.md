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
resolved as listed here. Implemented on `eftcamb-audit`: see "Status (implementation)" below.

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

## Status (implementation, 2026-10-21)

All approved items are committed on `eftcamb-audit` (`git log 3a2a4bb..HEAD`), never pushed. Final validation:
every agent-group boundary and all A9 commits build (15 commits checked in clean copies, plus HEAD);
`camb.tests.camb_test` passes (22 tests); pre-commit passes on all changed files; new `camb/tests/eftcamb_test.py`
(10 tests, ~11 s) passes; full regression and the result of every reproducer on HEAD are in `REGRESSION.md`
("repro" below = reproducer result on HEAD; "fixed" means the expected behaviour is observed).

| Item | Commits | Verification on HEAD |
|---|---|---|
| A1 | e20b307, 0e66c94 | repro eft-par-cache-realloc-abort/-crash fixed (get_BAO, reused Params, never-on RGR turn-on accepted) |
| A2 | 936a0a2 | repro eft-pickle-repr-loses-model, pickle-deepcopy fixed; eftcamb_test |
| A3 | cce814e, bde1b86 | getters/get_eft_functions/transfer-only evolution: values or CAMBValueError, no segfault |
| A4 | d38102a, d629225, 40444c2 | 3+1 background solver works; audit inis refused by the physics checks (open question 1) |
| A5 | e2a442b | repro fixed; eftcamb_test |
| A6, A22 | 0d36074 | repros fixed (model kept, EFT-only update); eftcamb_test |
| A7 | 785ce47 | repro fixed (47 keys, identical re-read); eftcamb_test |
| A8 | bc0324e | repro fixed |
| A9 / C12 | 3c759e7, 073de3d, 3a787c0, cf87fca, e756367, aa9a9a7, d7a9cda, 4819c8f | each commit builds; repros stop-kills-python, initialize-parameters-bad-value fixed (CAMBError/CAMBValueError with the Fortran message); all 85 runs bit-identical to F2 |
| A10 | aaf358d | checked build passes the stability check of 3_RPH_KT_power_1 |
| A11 | 1463967 | checked command-line exe of HEAD runs 1_EFT_GR past SetParams |
| A12 | 13b4f3d, cce9832 | repros fixed; 5_hdsk_TG_2 V_scf 23 (garbage) -> 0 |
| A13 | 7f96bad | repro fixed (ratio 1) |
| A14 | d983709, 42bbce3, 6fa6a74 | IC 3 = IC 1 for Omega = Gamma3 = 0; IC 2 = IC 1 |
| A15 | 878b118, df6840b, 6fe65d3 | turn-on 0, <0, NaN, >1 raise CAMBValueError, also with skip_RGR |
| A16 | 6ce961f, 0113598, 5b214ea | l_max_scalar 2000 with turn-on 1e-10 runs; 5b214ea: turn-on <= ~1e-24 segfaulted (found by this validation), now refused below tau = 1e-12 Mpc |
| A17 | 9e7f38c, 51b0327 | repro fixed (H0 to 4e-6) |
| A18 | 7a9bd1d | repro fixed |
| A19 | 1eb9367 | repros fixed (bit-identical after evolve_hubble = T) |
| A20 | 0de6891 | tensor-only run clean with FP traps on the checked build |
| A21 | 7a9bd1d, d0a1232, c7e004c, 3c41651 | string getters round-trip (eftcamb_test, read-parameters repro) |
| A23 | ca7de92, 14473b6, e024587, ccf65e1, edddf4f, 5fe97db, 6fa6a74 | repros rph-param-index, omp-private, qs-mu-nan fixed |
| B1 | 4209b57 | code verified; K-mimic single-mode glitch gone |
| B2 | 2d7e5ef | repro fixed (code/analytic = 1, ghost rejected) |
| B3 | a072620 | output-pidotdot repro fixed; Horava SolSyst_3 TT 0.40 at this commit |
| B4 | a78272a | both sampler repros fixed |
| B5 | 5fe97db | QS_sigma = general formula |
| B6, B12 | 27404a1, a26acac, 2d63421 | Hdotdotdot and additional derivatives vs finite differences 1e-4-7e-4; Minkowski verdict consistent |
| B7 | 28228ea | second Friedmann residual 1e-16; Horava spectra moved (REGRESSION.md) |
| B8, B15, B23 | d265ee6, 1423196, aba34f6, ec7bd8c, e5bed71 | Scaling Cubic H = CAMB H to 4e-16; Extended Galileon runs (was error), tracker background exact |
| B9 | f555f07 | EFTflag = 0, omk = +-0.05 vs CAMB 3e-14 |
| B10 | 9dc2900 | counts C_l vs CAMB 1.4e-7 |
| B11 | 97d774f, fe411a1 | no "dtauda is NaN"; OL-gamma GR point now refused as ghost (open question 6) |
| B13 | a37145e, 6493412, 90016be, 052eb5d, 7203f81, b1e1cf5, 748fa8e | all parametrization repros exact; numerical w_DE integral tabulated (runtime back to EFTwDE = 1 level) |
| B14 | 5db2def, efda336 | long keys read |
| B16 | 1423196, 3b63d3a, 9f944c7 | Beyond Horndeski H(z=0) = H0 (was 0.52 H0), spectra O(1) changes |
| B17 | 978cd52 | ini thetastar solved to 6e-8 for EFT models |
| B18 | b3e8280 | H from densities = H to 2e-16 |
| B19 | 8b11691 | extrapolation identical to upstream |
| B20 | 4f1decb | theta(T)/theta(F) - 1 = -2.5e-8 |
| B21 | 38eea73, 643f762, c6b0fb0 | 10x tables default (K-mimic included); H(z) changes 2-8e-5 |
| B22 | 42bbce3, 6fa6a74 | IC 2/3 series only where valid; IC 2 = IC 1 for late turn-on |
| C1 | 562e026 | omk != 0 with EFT: CAMBError |
| C2 | f152a94, 3ecaae6, 15370af, 4ebb732 | IC poles gone (RPH_MB_const_1, K-mimic); K-mimic turn-on dependence 1e8 -> 4% at k = 1 |
| C3 | a87186a, df6840b, 6fe65d3; tuned in da3c2d1 (Q3) | RGR threshold convergence; 4_Des_F_R_LCDM_2 P(k) -13% at k >= 1, runtime 40 -> 110 s; after Q3 (probe k = 1/Mpc, threshold 5e-4): turn-on 0.039, 36 s vs A0 33 s, P(k <= 1/Mpc) 4e-4 from the converged result |
| C4 (b)+(c) | 58cd024 | mass check prints "experimental" and is ignored; removed from the 14 inis; (a) left to the physics workflow |
| C5 | 2061238 | evolve_metric_h = T raises |
| C6 | d28a30a, ee26ecc, 2d6116d | 21cm refused; Halofit/HMcode get the background CPL (w0, wa): Pure EFT wCDM vs fluid P_nl/P_lin 3.3e-2 -> 5.5e-6 |
| C7 | e5d2b33 | Horndeski H0 adopted (5_hdsk_TG_1/2 sigma8 changes by 14%, now from the EFT H0) |
| C8 | cdcc3e3, 40444c2, 3c41651 | healthy RPH accepted at stability time 1e-10; threshold 1e-8 in the shipped inis |
| C9 | 2e86262 | turn-on jump sampled; low-l TT windows sensitivity 4.6e-2 -> 1.4e-4; RK45 step limit; get_evolution tolerance (partial, open question 4) |
| C10 | c9e1499 | README EFTCAMB section and MOD comments; EFTflag = 0 vs CAMB = roundoff with max_redshift = 50 |
| C11 | 006e9cf | set_dark_energy / set_classes raise with an active model |
| D | 209fdb5, f37f5f3, cce9832, d629225, d7bac4f (F-124), e5bed71 (F-106), b1e1cf5/90016be (F-116, F-117, F-109), ccf65e1 (F-096), 51b0327 (F-139), a87186a/df6840b (F-122) | done as listed; not done (S3, no approved change): F-108 exact real comparisons, F-121/F-131 massive-nu pressure approximation, F-123 test inis that run LCDM, F-126/F-127 sampler coverage, F-130/F-138 dead code, F-137 Makefile flags; F-135 documented only (C10); upstream-only F-099, F-100, F-115, F-120 untouched |

### Open questions for the user

1. **3+1 designer (007p8 audit inis).** A4 fixes the background solver, d629225/40444c2 accept the exact GR point,
   but both audit inis are still refused ("theory unstable": gradient instability -1.9e-26 at a = 1e-10 with threshold
   1e-8, beyond-GLPV warning with the ini's threshold 0). Is that physical? Also: when b3 = rho2/sigma1 -> 0 with
   rho2 != 0 somewhere, should the zz ghost test use the limit of Azz instead of skipping it?
2. **K-mimic turn-on sensitivity.** After C2, sigma8 still moves by ~2e-4 (0.61604 vs 0.61609-0.61625) and P(k) by
   1.5-3e-2 at k ~ 10/Mpc with the turn-on time or table density; the C2 branch chosen at the turn-on differs.
   Acceptable, or should the physics workflow look at K-mimic ICs?
3. **4_Des_F_R_LCDM_2 runtime** 40 s -> ~110 s (C3: earlier, converged turn-on 0.016 where the f(R) pi equation is
   stiff). Accept, or tune EFT_RGR_qs_threshold (1e-3) for speed?
4. **C9 residuals.** (a) GR-limit consistency: a near-GR EFT model (Omega = 1e-9, turn-on 0.01) now differs from
   EFTflag = 0 by TT 2.3e-5 (A0: 2.7e-6), from integrating exactly to the turn-on and the extra source samples;
   this is the ~2e-5 floor seen in most models. (b) 2_PEFT_Omega_const_1 low-l TT still changes by 2.1e-3 with
   TimeStepBoost = 3 (was 1.9e-2; GR 6e-5). (c) get_evolution at default accuracy: pi'(k ~ 1/Mpc) still O(1) off the
   IntTolBoost = 8 result (pi 3e-2). Accept, or refine further?
5. **Beyond Horndeski.** No massive neutrinos in its background (not approved, still missing); its hpprime contains
   an `OmegaRadprime` that looks like a typo (to be checked in the physics workflow).
6. **OL-gamma GR limit.** With B11 the "dtauda is NaN" rejection is gone, but 007p5 with OLOmega0 <= 1e-6 and its exact
   GR point are refused as ghosts at a = 1e-10 (kinetic -5e-3 / -9.4e-3 relative): the designer background leaves
   c != 0 at early times (c ~ -2.6e-4 at a = 1e-9). The shipped 007p5 inis (Omega0 = 0.05) run.
7. **3_RPH_MB_const_1 IC regularisation tolerance.** The C2 weight (0.1 of the term scale) removes the pole; P/P_GR
   still steps from 1.37 to 0.29 across k_c = 1.6e-3/Mpc (default and boost2 agree). Is the branch choice right?
8. **Upstream CAMB bug:** `get_background` (and Thermo_Init) with source windows and WantCls = False segfaults
   (bounds error on State%Redshift_W in the checked build); identical in plain CAMB v2 (B0). Report upstream?
9. **read_ini(no_validate=True)** with an unreadable EFTCAMB value still stops the interpreter (as for CAMB keys,
   upstream behaviour); the validated read_ini raises CAMBValueError. Change the CAMB ini reader too?
10. **F-135** (transfer-time mismatch error only with print_fortran_warnings) is documented only (C10), not changed.
11. **Horndeski theta.** With C7 (EFT H0 adopted, Horndeski_shooting = F), cosmomc_theta does not depend on the input
    H0, so theta -> H0 solvers cannot work for such models: refuse theta input for them?
12. **EFTCAMB_parameter_cache on input parameters** holds uninitialised values until the first calculation (it is
    filled in SetParams; harmless, but visible from Python). Initialise it at allocation?
13. Turn-on times below tau = 1e-12 Mpc (a ~ 2e-18) are now refused (5b214ea); previously a <~ 1e-22 ran and
    a <~ 1e-24 crashed. Fine as a hard limit?
14. The physics workflow (track 7, including the C4 (a) mass-eigenvalue re-derivation) has not run.

### Decisions on the open questions (user, 2026-10-21)

| # | Decision | Action |
|---|---|---|
| 1 | Leave as is | Flag for later: check the 3+1 designer stability conditions against the literature |
| 2 | Accepted | none |
| 3 | Tune for performance; accept up to ~20% slower than A0 (linear physics; small-scale P(k) less important) | **Done, da3c2d1:** quasi-static probe k = 10 -> 1/Mpc, `EFT_RGR_qs_threshold` 1e-3 -> 5e-4. 4_Des_F_R_LCDM_2 turn-on 0.016 -> 0.039, runtime 94 -> 36 s (A0 33 s, 1.08x; 1 thread, run alone). Vs an early turn-on (0.007): P(k <= 1/Mpc) 4.0e-4, P(k <= 10/Mpc) 3.3e-2, sigma8 5.4e-6, lensed C_l 2.8e-5, C_l^phiphi 1.1e-4 (A0: 1.6e-2, 0.15, 2.3e-4, 3.3e-5, 3.7e-3). No other test model or audit ini changes its turn-on; trade-off table in REGRESSION.md |
| 4 | Accepted | none |
| 5 | Leave as is | Flag for later: full Beyond Horndeski model check from the literature (massive neutrinos, `OmegaRadprime` term) |
| 6 | Accepted | none |
| 7 | Option (a): keep the C2 regularisation tolerance 0.1 | none; the physics workflow's re-derived mass check (C4 a) decides whether such tachyonic-at-turn-on models are rejected |
| 8 | No upstream report | none |
| 9 | No upstream CAMB changes | none |
| 10 | Accepted | none |
| 11 | Yes | **Done, eff41e9:** both theta -> H0 solvers (`CAMB_SetH0ForTheta`, ini path; `set_H0_for_theta`, Python `set_cosmology`) fail when the computed H0 differs from the trial H0: "thetastar/cosmomc_theta cannot be used with this model: its H0 (85.5373) is set by the EFT background" (CAMBValueError; CAMBError for `read_ini(no_validate=True)`). Horndeski_shooting = T also refused (the background still fixes H0); other EFT models (and designer-potential Horndeski, which starts from the input H0) unchanged; test `test_theta_input` |
| 12 | Yes | **Done, af2504f:** default initialisation (0, arrays unallocated) of every `TEFTCAMB_parameter_cache` component in the type, so every allocation (CAMB_SetDefParams, read_ini, command line, Python `F2003Class_new` with mold) is defined; test `test_parameter_cache_initialised` (CAMBparams(), new instance, read_ini: all 0); spectra bit-identical |
| 13 | Option (a): keep the hard error for turn-on earlier than tau = 1e-12 Mpc | none |
| 14 | Accepted (physics workflow later) | none |

Items 3, 11, 12 implemented on 2026-10-21 at the user's request (before the commit review): da3c2d1 (Q3), af2504f (Q12),
eff41e9 (Q11). Regression of all 85 inis vs F2: only 4_Des_F_R_LCDM_2 changes (Q3); `camb_test` + `eftcamb_test` pass.
Agent worktrees are kept.

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
