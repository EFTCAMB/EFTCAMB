export const meta = {
  name: 'eftcamb-audit-tracks',
  description: 'EFTCAMB deep audit: 15 track agents (layout, crash hunting, determinism, static sweep, GR limits, parametrizations, derivatives, robustness, stability, CAMB v2 integration, wrapper, regression, convergence) with adversarial verification of every finding',
  phases: [
    { title: 'Audit', detail: 'one agent per audit track, read-only on frozen builds' },
    { title: 'Verify', detail: 'adversarial reproduce/refute per finding' },
  ],
}

const AUDIT_ROOT = '/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit'

const CTX = `You are one agent in a multi-agent deep audit of EFTCAMB (effective-field-theory dark energy / modified gravity patch of the CAMB Boltzmann code; Fortran 2003 + Python ctypes wrapper). Your final answer is the structured result only.

REPO: /Users/marco/Desktop/IW/eftcamb_update/EFTCAMB, git branch eftcamb-audit at commit 5bef49a = the code under audit. EFTCAMB is one squashed commit "Applied EFTCAMB modification to CAMB" on top of upstream CAMB v2 (branch \`camb\`, 1d1d3d6). The pre-rebase EFTCAMB (on CAMB 1.6.5, upstream base 1db4141) is branch eftcamb-pre-rebase (76b5698). EFTCAMB sources: fortran/eftcamb/ (01 defs, 02 utilities/algorithms, 03 cache, 04/05 parametrizations, 06 abstract models, 07f designer models, 08f full-mapping models, 09 main/IC/RGR/stability, 10 background output, 11 python bridge). Hooks into CAMB are marked "EFTCAMB MOD START/END" in fortran/*.f90, the Makefiles and camb/*.py. Python wrapper: camb/eftcamb.py <-> fortran/eftcamb/11_EFTCAMB_python.f90. Test inis: fortran/eftcamb_test/parameters/*.ini (73 models + base_params.ini, hdsk_base_params.ini).

HARD RULES
- Do NOT edit tracked files in the repo, do NOT build in the repo checkout, do NOT run git commit/checkout/stash/reset, never push. You may create files only under fortran/eftcamb_test/audit/ (scripts, inis) and under ${AUDIT_ROOT}/work/<your-track>/ (outputs, tree copies).
- Scripts are Python (or Fortran test programs); no bash scripts. Put them in fortran/eftcamb_test/audit/checks/<your-track>/. Minimal reproducer for each S1/S2 finding: fortran/eftcamb_test/audit/repro/<finding-slug>.py, standalone, runnable with /Users/marco/Software/miniforge3/bin/python, printing observed vs expected numbers. Run it before reporting.
- Frozen builds under ${AUDIT_ROOT}/trees/ (never modify or rebuild them):
   A0  = code under audit, release flags (-O3)
   AD  = same code, gfortran -g -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-999999 -finit-logical=true -fbacktrace (no -DDEBUG)
   ADD = AD flags plus -DDEBUG (DebugEFTCAMB=.true.)
   B0  = plain upstream CAMB v2 (branch camb), release
   B1  = pre-rebase EFTCAMB on CAMB 1.6.5, release
  Each is a full source tree with camb/camblib.so; each tree's AUDIT_TREE.json gives ref/commit/flags. Build dirs (modules, objects) are ${AUDIT_ROOT}/build/<tree>/ (for A0: build/A0/camb/Releaselib and build/A0/forutils/...; look).
- For a modified build (instrumentation, candidate fix, patched debug build): cp -R ${AUDIT_ROOT}/trees/<T> ${AUDIT_ROOT}/work/<track>/<name>, edit the copy, then from the copy run: CAMB_BUILD_DIR=${AUDIT_ROOT}/work/<track>/<name>_build /Users/marco/Software/miniforge3/bin/python setup.py make (release flags; ~2-4 min). For non-default flags: cd <copy>/fortran && CAMB_BUILD_DIR=... make python PYCAMB_OUTPUT_DIR=../camb/ FFLAGS="<flags>" (see setup_trees.py for the flag strings; COMMON = -MMD -cpp -ffree-line-length-none -fallow-argument-mismatch -fmax-errors=4 -fopenmp -lblas -llapack). harness.run_model accepts a tree path.
- Shared machine: 18 cores, ~15 agents running at once. Use OMP_NUM_THREADS<=2 and at most 4 concurrent CAMB processes unless your check needs otherwise. A Bash call times out after 10 minutes: split long batches into chunks.
- Python: /Users/marco/Software/miniforge3/bin/python (3.13, numpy, scipy, sympy). gfortran 16.1 (conda-forge). macOS arm64.

HARNESS fortran/eftcamb_test/audit/checks/harness.py (import: sys.path.insert(0, "<repo>/fortran/eftcamb_test/audit/checks"); import harness):
  run_model(tree, ini, out_npz, overrides={ini_key: value}, evol_k=[k,...], threads=2, timeout=1800, python_set={"EFTCAMB.attr": v}, extra_env={}) -> status dict with status "ok" | "signal:SIGFPE" | "error:rc=N" | "timeout", log path, log_tail. Runs in a subprocess (cwd = output dir) importing camb only from the chosen tree; wraps the ini as DEFAULT(model.ini) + overrides (overrides win); forces EFTCAMB_write_background=F.
  run_batch(tree, inis, outdir, jobs=6, suffix="", **run_model_kwargs); model_inis(pattern="*.ini"); load(npz) -> dict; compare(a, b, lmax=2500) -> {cl:{spectrum:{TT,EE,BB,TE}}, pk:[per z], sigma8, bg_H, bg_tau, eft:{field: max|diff|/max|field|}, worst_cl, worst_pk, bit_identical} (TE normalised by sqrt(TT EE)); parse_ini(path) -> flattened dict.
  npz keys: cl_total, cl_unlensed_scalar, cl_unlensed_total, cl_lensed_scalar, cl_tensor, cl_lens_potential (raw muK^2, rows = ell, cols TT EE BB TE), pk_k [1/Mpc], pk_z, pk, sigma8, bg_z, bg_H, bg_tau, bg_DA; eft_a_grid (a = logspace(-8,0,401)) and eft_<field> for every EFTCAMB_timestep_cache field (get_eft_functions); evol_k, evol_a (logspace(-6,0,300)), evol_<field>[k, a] (get_scale_evolution) when evol_k given; meta (json: eft_params, eft_model_name, EFTflag, derived, errors, timings).
  CLI: python harness.py run|batch|compare --help.
  Reference runs at default settings, all 73 models: ${AUDIT_ROOT}/runs/A0/<model>.npz and ${AUDIT_ROOT}/runs/B1/<model>.npz (+ .json status, .log). (If a B1 file is missing, the batch may still be finishing; check again or rerun that model.)

KNOWN FACTS (already filed; do not re-report, build on them)
- F-001: Python calls the character(len=200) getters in 11_EFTCAMB_python.f90 (get_model_name, get_read_parameters, param names/labels, set_model_name) without gfortran's hidden string-length argument; -fcheck rejects it (11_EFTCAMB_python.f90:248, "actual string length 127/200"). The harness skips those calls on AD/ADD.
- F-002 (upstream CAMB): fortran/results.f90:2168 \`if (allocated(this%tb) .and. this%nthermo /= size(this%tb))\` evaluates size() of an unallocated array (Fortran .and. does not short-circuit) and aborts every -fcheck run in thermo init; same line exists in branch camb. Any -fcheck run needs a tree copy with this line rewritten as a nested if.
- B1's camb.read_ini ignores EFTCAMB keys; the harness then initialises EFTCAMB from the flattened ini (meta.eft_init_fallback=True).
- All 73 models run without error on A0 at default settings.
- EFTflag=0 A0 vs B0 differ up to 2.4e-5 (EE); suspected cause: EFTCAMB changed reionization max_redshift (CAMB default 50). Track 3 proves/disproves.
- CAMB v2 changed tensor Cls by ~14% (and scalar Cls by ~0.2-1%) relative to B1 even in GR (EFTflag=0): differences of that shape are inherited from CAMB, not EFT.
- Known not converged at default accuracy in both versions: 3_RPH_MB_const_1, 5_Kmimic_2, 5_Horava_SolSyst_*.
- Models without a test ini: Omega_Lambda_gamma (007p5), Omega_Lambda_alpha (007p6), ShiftSym_alphaB (007p7), 3+1 wDE designer (007p8), Beyond Horndeski (008p5), Extended Galileon (008p7). If your track needs them, create inis in fortran/eftcamb_test/audit/inis/ (find keys in the model's read_model_selection / init_model_parameters_from_file; DEFAULT(<abs path>/fortran/eftcamb_test/parameters/base_params.ini)). Reuse inis another agent already created there if present.

FINDINGS
- A finding is a defect in the code under audit (EFTCAMB code, or EFTCAMB's modifications of CAMB). Upstream-CAMB-only defects: put "upstream:" in the title, severity S3 unless they corrupt EFT results.
- Severity: S1 = wrong physics in a supported model, or crash/hang/memory corruption reachable with valid input; S2 = wrong only in an edge case/parameter region, or silently poor accuracy; S3 = robustness, hygiene, docs, dead code.
- Every S1/S2 needs a repro you actually ran; quote its numbers in observed. Group repetitive hygiene items (same pattern) into one finding listing all locations.
- Do NOT fix the repo. Propose a fix (unified diff or precise description). You may test it in your own tree copy and give before/after numbers. fix_class: "safe-local" (crash/bounds/uninit/leak/race, output change <= 1e-5), "physics-changing", or "unknown".
- Be skeptical of your own results before reporting: rule out script bugs, unit/normalisation mismatches, accuracy noise, known facts, intended documented behaviour.
- In checks[], record every check you ran with the measured value and threshold, PASS or FAIL. In not_covered[], list anything in your brief you could not do and why. Your verdict is PASS only if every check in your brief passed and nothing is uncovered.`

const FINDINGS_SCHEMA = {
  type: 'object',
  properties: {
    track: { type: 'string' },
    verdict: { type: 'string', enum: ['PASS', 'FINDINGS', 'INCOMPLETE'] },
    summary: { type: 'string' },
    commands: { type: 'array', items: { type: 'string' } },
    checks: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          name: { type: 'string' },
          value: { type: 'string' },
          threshold: { type: 'string' },
          pass: { type: 'boolean' },
          notes: { type: 'string' },
        },
        required: ['name', 'value', 'pass'],
      },
    },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          slug: { type: 'string' },
          title: { type: 'string' },
          severity: { type: 'string', enum: ['S1', 'S2', 'S3'] },
          file: { type: 'string' },
          line: { type: 'integer' },
          other_locations: { type: 'array', items: { type: 'string' } },
          repro_script: { type: 'string' },
          repro_command: { type: 'string' },
          observed: { type: 'string' },
          expected: { type: 'string' },
          root_cause: { type: 'string' },
          proposed_fix: { type: 'string' },
          fix_class: { type: 'string', enum: ['safe-local', 'physics-changing', 'unknown'] },
          confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['slug', 'title', 'severity', 'file', 'observed', 'expected', 'root_cause', 'proposed_fix', 'fix_class', 'confidence'],
      },
    },
    not_covered: { type: 'array', items: { type: 'string' } },
    artifacts: { type: 'array', items: { type: 'string' } },
  },
  required: ['track', 'verdict', 'summary', 'commands', 'checks', 'findings', 'not_covered'],
}

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    reproduced: { type: 'boolean' },
    refuted: { type: 'boolean' },
    verdict: { type: 'string', enum: ['CONFIRMED', 'PLAUSIBLE', 'REFUTED'] },
    corrected_severity: { type: 'string', enum: ['S1', 'S2', 'S3', 'not-a-bug'] },
    root_cause_agrees: { type: 'boolean' },
    reason: { type: 'string' },
    evidence: { type: 'string' },
  },
  required: ['reproduced', 'refuted', 'verdict', 'corrected_severity', 'root_cause_agrees', 'reason'],
}

const TRACKS = [
  { key: 't0-layout', prompt: `TRACK 0 (ctypes layout). Every later check reads Fortran data through ctypes, so verify the layouts exactly.
For EFTCAMB_timestep_cache and EFTCAMB_parameter_cache (Fortran in fortran/eftcamb/03_abstract_EFTCAMB_cache.f90; Python in camb/eftcamb.py), TEFTCAMB (Python class EFTCAMB in camb/eftcamb.py vs type TEFTCAMB in fortran/eftcamb/09_EFTCAMB_main.f90) and the EFT-related part/tail of CAMBparams (camb/model.py _fields_ vs type CAMBparams in fortran/model.f90): write a Fortran program compiled against the A0 build's module files and linked against the A0 objects or camblib.so that prints storage_size of each type and the byte offset of every component (loc(x%comp) - loc(x)); generate it from the Python field lists so nothing is skipped. Compare with ctypes sizeof/offsetof from the classes (import camb from the A0 tree). Check component order, types and kinds (Fortran default logical is 4 bytes vs c_bool 1 byte: find how CAMB's baseconfig maps logicals and whether EFTCAMB follows it; integer kinds; allocatable array descriptors via AllocatableArrayDouble). Then a value round trip: for 3 models (designer f(R), Horndeski hdsk_TG, K-mouflage) compare a few fields returned by get_eft_functions with independent values (adotoa vs a*H from results.hubble_parameter in CAMB units; grhoc_t, grhob_t vs the analytic a-scaling of the parameter-cache densities; EFTOmegaV at a=1 for f(R) vs the B0 relation; a field echo). Check get_evolution/get_scale_evolution: shape/order of the returned array ([k][time] vs Fortran fill order in EFTCAMB_GetEvolution in 11_EFTCAMB_python.f90), and that np.sort(times) inside get_evolution cannot misalign results from the caller's eta order (it sorts silently). Any mismatch that corrupts values read by Python is S1.` },
  { key: 't1-crash', prompt: `TRACK 1 (crash hunting with runtime checks). Find every runtime error, floating-point trap, bounds violation, use of uninitialised data, NaN/Inf in EFTCAMB.
(1) Copy tree AD to your work dir, patch the upstream results.f90 size() line (F-002) into a nested if, rebuild with the AD flags (from trees/AD/AUDIT_TREE.json). If further upstream-only trips appear, patch each in your copy and record them together as one grouped "upstream:" S3 finding with all locations; keep going until EFT code is exercised.
(2) Run all 73 test models through Python on that build (harness.run_batch with your tree path, threads=1, jobs<=6, timeout large; -O0 is slow, chunk the batch). For every non-ok status read the log, get the backtrace (file:line), decide EFT vs upstream, minimise to a repro. Also scan outputs of ok runs for NaN/Inf.
(3) Build the command-line executable in your copy (cd fortran && make camb with the same FFLAGS and CAMB_BUILD_DIR) and run it from a scratch cwd on at least one model per family (PEFT, RPH, Des_F_R, Des_mc5e, ADE, hdsk_TG, Horava, Kmimic, Kmouflage, quint_*, scg) plus every model that failed in Python; note differences between the two entry points.
(4) Make a patched copy of ADD the same way and run a handful of models to exercise the DebugEFTCAMB path (cmbmain.f90 opens cache files on the master EV but the OMP source loop uses PRIVATE(EV); there is a \`stop\` if fixq_array is set). Report defects in that path.
(5) Run the 6 models that have no test ini (create minimal inis, see KNOWN FACTS) on the checked build.
(6) Compiler diagnostics: compile fortran/eftcamb sources (and the EFT MOD blocks' files) with -Wall -Wextra -Wconversion -Wimplicit-interface -Wuninitialized -Wmaybe-uninitialized -O2 in your copy (note the release Makefile passes -w which hides all warnings); triage: list individually warnings that indicate real bugs (uninitialised use, precision-losing conversions in physics, character truncation, unused dummy that should be used, unreachable code), group the rest by category with counts. Compile fortran/eftcamb/08f_full_models/hdsk_coefficients*.c with -Wall -Wextra and check the Fortran interfaces that call them (bind(c) names, argument types/kinds, value vs reference, array sizes).` },
  { key: 't2-state', prompt: `TRACK 2 a-c (determinism and state). (The static sweep 2d is another agent.)
(a) Thread determinism: for at least 30 models covering every family, run A0 with OMP_NUM_THREADS=1 and =8 (8-thread runs at most 2 at a time) and compare (harness.compare + exact comparison of eft_* arrays). Establish the baseline first: same comparison for 1_EFT_GR and on B0. Any relative diff > 1e-10 beyond the GR baseline is a race candidate: repeat same-thread-count runs to separate run-to-run nondeterminism from thread-count dependence; localise (EFT functions vs perturbations; which k; which variables in get_scale_evolution).
(b) State leakage: in ONE Python process importing camb from A0, compute model A, then B, then A again, for pairs: EFTflag 0 -> nonzero -> 0; designer -> full mapping; two full-mapping models; two parametrizations of the same model; with/without massive neutrinos; a model that fails (unstable / error) between two good runs. Use read_ini + get_results and also the set_params path; also reuse one CAMBparams object after changing EFT parameters, copy.deepcopy a CAMBparams carrying EFTCAMB, and pickling if supported. The repeated A must be bit-identical to the first and to a fresh process. Localise any difference (which state persists) and find it in the Fortran (module variables, SAVEd locals, allocations not reset, designer caches/interpolation tables, RGR/stability module state, global_error_flag).
(c) Memory: 300 repeated get_results in one process, alternating 2-3 EFT models (also with read_ini each time, and with set_params each time); record RSS every 10 iterations (resource.getrusage / psutil if present / ps); report growth per iteration vs a GR baseline run the same way. If it grows, find the leak (initialize_parameters re-initialises the EFTCAMB object in place: are model components / parametrizations deallocated?).` },
  { key: 't2d-static', prompt: `TRACK 2d (static sweep: thread safety and numeric hygiene) over fortran/eftcamb/**/*.f90 and every EFTCAMB MOD block in fortran/*.f90. CAMB evolves k-modes inside OpenMP parallel loops (cmbmain.f90: source/transfer loops), so EFT code executed per k (equations.f90 derivs/initial/output -> EFTCAMB model compute_* procedures, timestep cache, pi factors, stability called per k?) must not write shared state. Find and assess:
- locals initialised in their declaration inside procedures (implicit SAVE): list all; flag those reachable from per-k evolution (race) or persisting across runs (state leakage);
- module-level variables and components of shared objects (CP%EFTCAMB, CP%EFTCAMB%model, parametrization objects) written during perturbation evolution or source output;
- stop / error stop in library code reachable from Python (kills the interpreter) instead of the error handler / global_error_flag;
- real ==, /= comparisons with physical meaning; integer division in physics expressions; real literals without _dl kind inside double-precision physics expressions (count per file; list only those that change results above 1e-8);
- intent(out) dummies not set on all paths; function results not set on all paths; uninitialised locals;
- division by quantities that can vanish (1+Omega, adotoa, k, a, Hdot - adotoa^2, denominators of stability/QSA/pi-equation factors) without guard;
- copy-paste errors (an index/name that should differ, e.g. Gamma3P where Gamma4P intended) across the parametrization and model files and in the pi/Einstein factor code;
- OpenMP directives in EFT code: shared/private correctness.
Read enough context to decide each candidate; for thread races confirm reachability inside the parallel region (follow the call chain). Produce grouped findings with every location. Where a race is plausible, try to demonstrate it (e.g. instrumented copy or comparing 1 vs many threads for the affected model).` },
  { key: 't3-gr', prompt: `TRACK 3 (GR / LCDM limits).
(a) EFTflag=0 on A0 vs B0 (plain CAMB v2): GR ini plus variants (massive neutrinos, omk != 0, tensors, lensing, transfer at several z, w != -1 fluid dark energy with EFTflag=0, accuracy_boost=2). Must be bit-identical except for intentional documented changes. They differ by up to 2.4e-5 in EE now. Find every cause: list the EFTCAMB MOD blocks active when EFTflag=0 (git diff camb eftcamb-audit -- fortran camb; note the reionization max_redshift change), and prove causality by building a copy of A0 with each candidate reverted until results are bit-identical with B0 (or list what remains). Any EFTflag=0 difference not explained by an intentional, documented change is a finding.
(b) For each model family find its GR/LCDM-equivalent parameter point (Pure EFT with all functions zero and w=-1; designer f(R) with tiny B0; designer mc quintessence with w=-1 if allowed; RPH with zero functions; Horava at GR values; Horndeski/K-mouflage/quintessence/ADE/scaling cubic/extended galileon where the scalar decouples, if such a point exists; Omega_Lambda_gamma/alpha with zero functions; Beyond Horndeski; 3p1; ShiftSym). Run each and compare with B0 LCDM at the same cosmology: report max deviation in Cls/P(k)/sigma8/H in checks; flag > 1e-4. Scan the approach to GR (e.g. B0 or the coupling = 1e-2, 1e-4, 1e-6, 1e-8, 0) and check the deviation goes smoothly to the EFTflag=0 result, with no jump from the RGR logic (EFTCAMB_GR_threshold, 09_EFTCAMB_RGR.f90) or from switching between "GR" and "EFT" code paths in equations.f90.` },
  { key: 't4a-param', prompt: `TRACK 4a (1D and 2D parametrizations). Units under test: every concrete type in fortran/eftcamb/04f_parametrizations_1D/*.f90 and 05f_parametrizations_2D/*.f90, through the abstract interfaces in 04_abstract_parametrizations_1D.f90 and 05_abstract_parametrizations_2D.f90, plus the allocator 04p1_parametrizations_1D_allocator.f90.
Write a Fortran driver compiled against the A0 build modules and linked against the A0 objects/library that instantiates each parametrization with several parameter sets (typical and edge values) and evaluates value, first, second, third (fourth if present) derivatives and the integral function on a grid in a (x for 2D) including a=1e-8, a=1, and points around internal nodes / turning points / steps. First find the exact definition of each method in the abstract interface (derivatives with respect to a? the "integral" is EFTCAMB-specific, e.g. for w(a) it may be exp(-3 int (1+w) dln a) or similar: determine it and test against that definition). Compare derivatives with high-order centred finite differences / Richardson extrapolation and integrals with scipy quad (in Python, from values printed by the driver). Tolerances: 1e-6 relative for first derivatives, 1e-4 for higher; splines: continuity of value and derivatives at nodes. Also check parameter set/get (set_param_array, parameter_value, parameter_names count matches number of parameters), feedback, init_from_file keys, and that each allocator selection integer maps to the documented type. Note the duplicate file prefix 04p019 (power_law_sum and spline5): check the build lists both. Report the worst error per method per parametrization in checks.` },
  { key: 't4b-derivs', prompt: `TRACK 4b (time derivatives stored in the EFTCAMB timestep cache). Using get_eft_functions (and get_scale_evolution where needed) on A0 for every test model (reference runs have eft_<field> on a=logspace(-8,0,401); write your own script importing camb from the A0 tree for denser grids), check every derivative relation numerically. First establish conventions from fortran/eftcamb/03_abstract_EFTCAMB_cache.f90 and the model code: CAMB uses conformal time tau; adotoa = conformal Hubble; Hdot = d(adotoa)/dtau etc.; EFTOmegaP/PP/PPP/PPPP are derivatives with respect to a; EFTcdot, EFTLambdadot are tau derivatives; Gamma_iP... w.r.t. a; alpha*dot; gpresdotm_t, gpresdotdotm_t, grhonudot_tot, gpinudot_tot, gpinudotdot_tot. Use high-order finite differences in ln a on smooth stretches (avoid turn-on and RGR switch times). Tolerance 1e-6 for first derivatives and 1e-4 for higher, relative to the field's scale; report the worst per field per model in checks (summarise). Also verify relations that must hold for any model: Friedmann (adotoa^2 vs total grho a^2/3 in CAMB units including the EFT dark-energy density), the Hdot relation with total pressure, matter conservation relations for the *_t densities. Report problems at turn-on/RGR times separately (they enter initial conditions).` },
  { key: 't5-robust', prompt: `TRACK 5 (numerical robustness). For representative models of every family (all models if time allows), vary EFTCAMB_back_turn_on, EFTCAMB_pert_turn_on, EFTCAMB_GR_threshold, EFTCAMB_stability_time (find defaults in the code) by a factor 10 each way and EFT_IC_type over all allowed options (1, 2, 3), measuring the change in Cls (ell<=2500), P(k), sigma8 vs defaults. First establish the numerical noise floor with accuracy_boost=2 vs 1 for the same models. Changes above 1e-4 (and above the noise floor) need a cause: physical (model not GR-like before turn-on, documented) or a defect: discontinuity of EFT functions at turn-on, initial conditions not on the attractor (09_EFTCAMB_IC.f90: pi, pi' initial values), RGR switch zeroing functions abruptly, designer interpolation tables starting late. Inspect pi(a), pidot(a) with get_scale_evolution around the turn-on for transients/oscillations indicating bad ICs. Also check EFTCAMB_evolve_delta_phi and EFTCAMB_evolve_metric_h on/off for models that support them (spectra must agree to < 1e-4) and EFTCAMB_use_background for Horndeski.` },
  { key: 't6-stability', prompt: `TRACK 6 (stability logic: fortran/eftcamb/09_EFTCAMB_stability.f90, each model's additional_model_stability, positivity bounds / Minkowski limit code, camb/camb.py positivity allowlist). For each model with stability checks: take the analytic conditions (scalar no-ghost, e.g. alpha_K + 3/2 alpha_B^2 > 0 in Horndeski; tensor no-ghost M*^2 > 0; gradient c_s^2 > 0, c_T^2 > 0; mass conditions) from the EFTCAMB papers (arXiv:1312.5742, numerical notes 1405.3590, alpha basis 1404.3713; Frusciante & Perenon review 1907.03150). Scan 1-2 parameters across the stable/unstable boundary finely enough to locate it to ~1% and compare the code's verdict (error raised / message) with the boundary computed independently from the cached EFT functions (rerun with EFTCAMB_skip_stability=T to get them) and the analytic formulas. Check: (i) a clean Python exception with a reason (no crash, no hang, process usable after); (ii) time sampling: build a case with an instability only in a narrow interval and see if it is missed; (iii) the checked interval covers the integrated interval (pert turn-on to today); (iv) EFT_mass_stability_rate and EFTCAMB_stability_threshold semantics; (v) ghost_math vs physical ghost flags; (vi) positivity bounds and Minkowski limit for the allowlisted models; is the allowlist consistent with which models implement the needed higher derivatives? Focus on decision logic; the cached stability coefficients' formulas are also checked by the physics workflow.` },
  { key: 't8a-equations', prompt: `TRACK 8a (CAMB v2 integration in fortran/equations.f90). EFTCAMB was rebased from CAMB 1.6.5 (eftcamb-pre-rebase, upstream base 1db4141) onto CAMB v2 (camb). List all EFTCAMB MOD START/END blocks in fortran/equations.f90 at HEAD (about 64). For each block: compare with the same block in eftcamb-pre-rebase, and compare the surrounding upstream code between 1db4141 and camb (git diff 1db4141 camb -- fortran/equations.f90) to see whether v2 changed the context: renamed or re-meaning variables, moved logic, new terms (new species, new source terms, changed y-vector indexing, massive-neutrino approximation switches, tight coupling, radiation streaming, output() sources, initial conditions), such that the EFT block now (i) misses a term EFT needs (e.g. a new contribution to dgrho/dgq/dgpi that must be passed to EFT), (ii) uses a variable whose meaning changed, (iii) sits in the wrong branch, or (iv) is duplicated or dead. Pay special attention to EV% indices of EFT variables (pi, pidot or delta_phi, h'), the replaced Einstein equations (z, sigma, etak), ISW/lensing source terms (EFTISW, EFTLensing), tensor equations, massive neutrinos, the dark-energy perturbation interface, and the tight-coupling switch. Give a per-block verdict (ok / suspicious / bug) in checks. Each bug needs a repro (a model where the difference matters; compare against B1 and against GR consistency).` },
  { key: 't8b-integration', prompt: `TRACK 8b (CAMB v2 integration outside equations.f90). Same method as for equations.f90 (compare each EFTCAMB MOD block with the pre-rebase version and with how upstream changed around it between 1db4141 and camb) for: fortran/cmbmain.f90 (13 blocks), results.f90 (12), camb.f90 (10), lensing.f90 (6), halofit.f90 (3), model.f90 (2), massive_neutrinos.f90 (2), reionization.f90 (1), fortran/Makefile and Makefile_main, and the Python side (grep -n "EFTCAMB MOD" camb/*.py). In addition: when EFTflag != 0, is every place where v2 computes background quantities (dtauda, conformal time tables, sound horizon, H(z), distances, DarkEnergy BackgroundDensityAndPressure, growth, sigma8, z_star/theta solvers, derived params, thermal history) routed through the EFT background (and EFTCAMB_use_background)? Look for v2 code paths that bypass it. Check the camb/results.py EFT mod (zero-crossing extrapolation of the matter power interpolation) for correctness, and the camb/model.py set_H0_for_theta EFT mod. Explain the reionization max_redshift change and its effect. Per-block verdicts in checks.` },
  { key: 't8c-features', prompt: `TRACK 8c (new or changed CAMB v2 features with EFT on). From git log / git diff 1db4141..camb and docs (docs/, CAMB_README.rst, CHANGELOG if present), list user-facing features added or changed in v2 and test each with EFT on (models: designer f(R), Horndeski hdsk_TG, K-mouflage, quintessence, a PEFT gamma model): automatic lens_potential_accuracy / kmax (lens_potential_accuracy=None default): does it under-resolve MG models? compare with lens_potential_accuracy=8; source windows / number counts / get_transfer_functions with EFT (are EFT-modified potentials used?); non-linear: halofit and HMcode with EFT on (consistent, or refused explicitly? check do_nonlinear in the ini path and NonLinearModel in Python); dark-energy model interface (EFTCAMB forces the fluid class: what if the user also sets dark_energy_model, w, wa with EFTflag != 0: double counting?); several massive neutrinos and the hierarchy setter; omk != 0 (supported or refused?); set_H0_for_theta / cosmomc_theta / thetastar with EFT; get_background, get_age, get_zre_from_tau; 21cm / redshift sources; changed accuracy defaults. For each: verdict in checks; defects as findings with repro.` },
  { key: 't9-wrapper', prompt: `TRACK 9 (Python wrapper: camb/eftcamb.py <-> fortran/eftcamb/11_EFTCAMB_python.f90, plus EFT parts of camb/camb.py, camb/model.py, camb/_ini.py, camb/results.py).
(i) get_eft_functions, get_evolution, get_scale_evolution return what Fortran computes: spot-check several fields for 3 models against independent routes (EFTCAMB_write_background=T text output written by a run in a scratch dir; or the ADD tree debug cache files; or a Fortran driver), with scalar inputs, arrays of k, unsorted eta, a before the perturbation turn-on, a=1.
(ii) read_parameters / boolify / autoconvert round trip for every test ini: ini -> read_parameters dict -> initialize_parameters on a fresh CAMBparams -> results bit-identical to the ini path. Look for conversion pitfalls: ints vs reals, T/F, strings, Fortran D exponents (1d-3) that Python float() rejects.
(iii) camb.set_params(**dict) with EFT keys vs read_ini: identical results.
(iv) write_ini (camb/_ini.py, CAMBparams.write_ini) with EFT on: are EFT keys written? ini -> read -> write -> read -> results identical? (Project rule: new class parameters must be handled in write_ini.)
(v) Error handling: invalid model selection numbers, missing parameters, unstable model, NaN/negative parameters: Python exception with message, no crash, process usable afterwards with results identical to a fresh process.
(vi) F-001 (missing hidden string-length argument): enumerate every affected interface, check whether release builds can overflow/garble (names > 200 chars, non-padded buffers) and propose the concrete fix (bind(c) with explicit lengths, or passing lengths from ctypes). Treat F-001 as known: report only new aspects.
(vii) initialize_parameters receives ALL set_params keys (CAMB keys too) and writes them into the EFT ini reader: can a CAMB key collide with an EFT key or change EFT behaviour?` },
  { key: 't10a-regression', prompt: `TRACK 10a (regression vs pre-rebase B1 and legacy B2). Reference runs: ${AUDIT_ROOT}/runs/A0 and ${AUDIT_ROOT}/runs/B1 (73 models, defaults). B2 = fortran/eftcamb_test/eftcamb_legacy/{benchmark_results,spectra_results} (text outputs of legacy EFTCAMB: read its README and file formats).
1. compare(A0, B1) for every model. Classify each diff > 1e-4: inherited from CAMB v2 if it has the same size and ell/k shape as the GR diff (runs/A0/1_EFT_GR vs runs/B1/1_EFT_GR): compute per ell/k the double ratio (A0/B1)_model / (A0/B1)_GR - 1; if |.| < 1e-4 it is inherited. Otherwise EFT-specific -> finding candidate. To shrink inherited differences you may rerun both trees with matched settings (higher accuracy, fixed lens_potential_accuracy, reionization max_redshift matched).
2. Compare with B2 legacy spectra where available (units, columns; legacy used older CAMB: compare EFT/GR ratios rather than raw spectra).
3. For every EFT-specific deviation, localise it: EFT functions differ (compare eft_* fields A0 vs B1) vs perturbations vs outputs; find the responsible change (git diff eftcamb-pre-rebase eftcamb-audit -- fortran/eftcamb camb/eftcamb.py, and the MOD blocks); prove causality with a copy of A0 that reverts the change when feasible.
Put a per-model table (worst inherited-corrected deviation for TT/EE/TE/BB/PP/P(k)) in checks.` },
  { key: 't10b-convergence', prompt: `TRACK 10b (numerical convergence). For every test model run A0 at accuracy_boost=2 and 3 (reference default runs exist; keep l_accuracy_boost and lens_potential_accuracy consistent, consider l_sample_boost) and measure default vs boost3 and boost2 vs boost3 for Cls ell<=2500 (TE normalised), P(k), sigma8. Use GR (EFTflag=0) at the same settings as the noise baseline. Flag models whose default-accuracy error exceeds ~3x the GR level or 1e-3. For non-converged models (known: 3_RPH_MB_const_1, 5_Kmimic_2, 5_Horava_SolSyst_*; plus any you find) find the root cause: which accuracy knob fixes it (ODE tolerance / time steps, k sampling of sources, ell sampling, transfer kmax, EFT-specific: designer background tables, B0 shooting tolerance, RGR/turn-on time, interpolation) and inspect the evolution (get_scale_evolution) for stiffness or under-resolved oscillations of pi (large mass, tiny or huge c_s^2). Propose a fix (EFT-specific step control, better ICs, tolerances) and test it in a tree copy if feasible.` },
]

const RESOURCES = `

RESOURCE RULES (mandatory; an earlier run overloaded the user's machine):
- Every CAMB run goes through harness.run_model / run_batch (they hold machine-wide CPU slots, run at nice 10, and exit if you die). Use threads=1 (2 only when needed); never more than 8 threads for a single run, and only for a determinism check. If you run CAMB from your own script instead, wrap the computation in \`with harness.slot(n_threads):\` and set OMP_NUM_THREADS=n_threads before importing camb.
- Rebuild tree copies only with harness.build_tree(src_tree, dest_dir, fflags=None|"<flags>", make_args=None) (one build at a time on the machine, make -j4); reuse copies from the earlier attempt when they are still valid instead of rebuilding.
- Never use run_in_background for CAMB, builds or long loops: run foreground chunks under 9 minutes (orphaned background jobs survived the last session and kept loading the machine).
- Delete large outputs you no longer need (tree copies of finished experiments, bulky npz) at the end; keep what your findings cite.`

function makeLimiter(n) {
  let active = 0
  const queue = []
  const next = () => {
    if (active >= n || !queue.length) return
    active++
    const job = queue.shift()
    job.fn().then(job.res, job.rej).finally(() => { active--; next() })
  }
  return fn => new Promise((res, rej) => { queue.push({ fn, res, rej }); next() })
}
const limit = makeLimiter(8)
const lagent = (prompt, opts) => limit(() => agent(prompt, opts))

function verifyPrompt(track, f, lens) {
  const lensText = lens === 'reproduce'
    ? 'LENS: reproduce and attribute. Run the repro script exactly as given, then write your own independent check of the same claim (different code path or method). Is the defect in the code under audit, upstream CAMB, or the test itself? Is the stated root cause (file:line) right?'
    : 'LENS: physics and expectation. Is the expected behaviour actually correct? Check against the equations/papers (EFTCAMB numerical notes arXiv:1405.3590, 1312.5742; alpha basis 1404.3713) and the code\'s documented intent. Could the observed behaviour be intended, a convention difference, accuracy noise, or a units/normalisation mismatch? Is the severity right?'
  return `${CTX}

You are an ADVERSARIAL VERIFIER for one finding reported by the ${track} agent. Try hard to refute it; confirm only what you can reproduce and cannot explain away. Do not fix anything. Write any scripts under fortran/eftcamb_test/audit/checks/verify/ and outputs under ${AUDIT_ROOT}/work/verify/.
${lensText}

FINDING:
${JSON.stringify(f, null, 2)}

Return: reproduced (did you see the effect yourself), refuted (is the claim wrong as a defect of the code under audit), verdict CONFIRMED (reproduced, real defect) / PLAUSIBLE (real-looking but not fully reproduced) / REFUTED, corrected_severity, root_cause_agrees, reason, evidence (numbers you measured).${track === 't0-layout' ? '' : RESOURCES}`
}

const results = await pipeline(
  TRACKS,
  t => (t.key === 't0-layout' ? agent : lagent)(`${CTX}\n\n${t.prompt}\n\nYour track id: ${t.key}. Work directory: ${AUDIT_ROOT}/work/${t.key}/.${t.key === 't0-layout' ? '' : `

RESUMING: earlier attempts of this exact track were interrupted (usage limits, one machine overload) after substantial progress; each attempt resumed from the previous one. Their condensed transcripts (notes, every command, short result tails) are ${AUDIT_ROOT}/prev/audit_${t.key}.md (attempt 1), ${AUDIT_ROOT}/prev/audit_${t.key}.run2.md (attempt 2) and, if it exists, ${AUDIT_ROOT}/prev/audit_${t.key}.run3.md (attempt 3). Scripts are in fortran/eftcamb_test/audit/checks/${t.key}/, outputs in ${AUDIT_ROOT}/work/${t.key}/, repro scripts in fortran/eftcamb_test/audit/repro/ (other tracks wrote some too). Read the most recent notes file fully first; consult earlier ones only as needed (grep them). Reuse sound results (rerun cheaply to confirm key numbers), do not redo finished work, and finish efficiently: aim to return your structured result well within your budget, listing anything left undone in not_covered rather than starting new long explorations. Repro scripts must still be run before you report them.${RESOURCES}`}`, {
    label: `audit:${t.key}`, phase: 'Audit', schema: FINDINGS_SCHEMA,
  }),
  (res, t) => {
    if (!res) return null
    const fs = res.findings || []
    log(`${t.key}: verdict ${res.verdict}, ${fs.length} findings (${fs.filter(f => f.severity !== 'S3').length} S1/S2)`)
    return parallel(fs.map(f => () => {
      // S3 findings are recorded unverified to save machine time and tokens
      if (f.severity === 'S3') return Promise.resolve({ ...f, track: t.key, verdicts: [], unverified: true })
      // one reproduce-lens verifier per S1/S2 finding, plus the physics lens for S1 (usage limits)
      const lenses = f.severity === 'S1' ? ['reproduce', 'physics'] : ['reproduce']
      return parallel(lenses.map(lens => () => lagent(verifyPrompt(t.key, f, lens), {
        label: `verify:${t.key}:${f.slug}:${lens}`, phase: 'Verify', schema: VERDICT_SCHEMA,
      }))).then(vs => ({ ...f, track: t.key, verdicts: vs.filter(Boolean) }))
    })).then(verified => ({ ...res, findings: verified }))
  },
)

return results.filter(Boolean)
