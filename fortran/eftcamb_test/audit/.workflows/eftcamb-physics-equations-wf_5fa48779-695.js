export const meta = {
  name: 'eftcamb-physics-equations',
  description: 'EFTCAMB track 7: independent derivations of the shared EFT equations (two derivations + reconciliation + Bianchi check), per-model derive-first checks for all 16 models, cross-code-path checks, with adversarial verification',
  phases: [
    { title: 'Prep', detail: 'inis and representative runs for all 16 models' },
    { title: 'Reference', detail: 'two independent derivations of the EFT equations, reconciled' },
    { title: 'Shared', detail: '7a checks on every model using the reference equations' },
    { title: 'Models', detail: '7b derive-first per model' },
    { title: 'Cross', detail: '7c same physics through different code paths' },
    { title: 'Verify', detail: 'adversarial verification of each finding' },
  ],
}

const AUDIT_ROOT = '/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit'
const PHYS = 'fortran/eftcamb_test/audit/physics'

const CTX = `You are one agent in a multi-agent physics validation of EFTCAMB (effective-field-theory dark energy / modified gravity patch of the CAMB Boltzmann code; Fortran 2003 + Python ctypes wrapper). Your final answer is the structured result only.

REPO: /Users/marco/Desktop/IW/eftcamb_update/EFTCAMB, git branch eftcamb-audit at commit 5bef49a = the code under audit (EFTCAMB as one squashed commit on upstream CAMB v2, branch \`camb\`). EFTCAMB sources: fortran/eftcamb/ (03 cache with EFTCAMB_timestep_cache fields; 04/05 parametrizations; 06_abstract_EFTCAMB_model.f90 with compute_Einstein_factors / compute_pi_factors / compute_tensor_factors / compute_stability_factors / compute_QS_mu / Sigma; 06p1 full mapping; 06p2 designer; 07f designer models; 08f full-mapping models; 09 main/IC/RGR/stability). EFT perturbation equations are used in fortran/equations.f90 (EFTCAMB MOD blocks). Python: camb/eftcamb.py (get_eft_functions, get_evolution, get_scale_evolution return arrays of EFTCAMB_timestep_cache). Test inis: fortran/eftcamb_test/parameters/*.ini.

HARD RULES
- Do NOT edit tracked repo files, do NOT build in the repo checkout, no git commit/checkout/stash/reset, never push. Create files only under fortran/eftcamb_test/audit/ and ${AUDIT_ROOT}/work/<your-id>/.
- Python or Fortran only (no bash scripts). Derivations: sympy scripts under ${PHYS}/<your-id>/ that run from scratch with /Users/marco/Software/miniforge3/bin/python and print their checks. Each check cites the paper and equation number it implements (arXiv id + eq.). Minimal reproducer per S1/S2 finding in fortran/eftcamb_test/audit/repro/<slug>.py (standalone; prints observed vs expected); run it before reporting.
- INDEPENDENCE: derive first from the action / source papers, save the derivation, and only then open the Fortran that computes the corresponding quantity. Never transcribe equations from the EFTCAMB code. You MAY read CAMB's GR code and the EFT code's call sites in equations.f90 to learn variable conventions and which equation a cached factor multiplies, but not the expressions that compute EFT quantities until your independent expression is written to disk (record in your result that you followed this order).
- Papers: fetch with WebFetch/WebSearch when available (arxiv.org/abs/<id>; TeX source at https://arxiv.org/e-print/<id> for exact equations). Core references: EFTCAMB implementation arXiv:1312.5742; EFTCAMB numerical notes arXiv:1405.3590 (defines the cached factor names EFTeom*, EFTpi*, tensor factors, in the version matching the code — check the latest arXiv version); alpha basis arXiv:1404.3713; QSA mu/Sigma arXiv:1606.05339; EFT review arXiv:1907.03150. Per model: papers cited in the model file. If a paper cannot be fetched, derive from first principles and state conventions explicitly.
- Frozen builds under ${AUDIT_ROOT}/trees/: A0 = code under audit (release); AD/ADD = debug variants (need a patch, avoid); B0 = plain CAMB v2; B1 = pre-rebase EFTCAMB. Never modify them. For a modified build, copy the tree into your work dir and rebuild there with CAMB_BUILD_DIR=<your work dir>/<name>_build python setup.py make (from the copy).
- Shared machine with ~30 agents: OMP_NUM_THREADS<=2, at most 3 concurrent CAMB processes. Bash calls time out after 10 minutes.
- Python: /Users/marco/Software/miniforge3/bin/python (3.13; numpy, scipy, sympy).

HARNESS fortran/eftcamb_test/audit/checks/harness.py (sys.path.insert(0, "<repo>/fortran/eftcamb_test/audit/checks"); import harness): run_model(tree, ini, out_npz, overrides={ini_key: value}, evol_k=[...], threads=2, python_set={...}) -> status; run_batch; model_inis(pattern); load(npz); compare(a, b) (Cls TT/EE/BB/TE with TE/sqrt(TT EE), pk per z, sigma8, eft fields; worst_cl, worst_pk, bit_identical); parse_ini. npz: cl_* (raw muK^2), pk_k [1/Mpc], pk_z, pk, sigma8, bg_z, bg_H, bg_tau, bg_DA, eft_a_grid (logspace(-8,0,401)) + eft_<field> (get_eft_functions), evol_k, evol_a (logspace(-6,0,300)) + evol_<field>[k,a] (get_scale_evolution) when evol_k given; meta json. Reference runs (defaults, 73 test models): ${AUDIT_ROOT}/runs/A0/<model>.npz. For finer grids or other outputs, write a script that imports camb from ${AUDIT_ROOT}/trees/A0 (sys.path.insert(0, tree)) and calls cp = camb.read_ini(ini); results = camb.get_results(cp); cp.EFTCAMB.get_eft_functions(results, a); cp.EFTCAMB.get_scale_evolution(results, k_list, a) (returns (fields, structured array [k][time])). Run such scripts in a subprocess or separate process per model.
- Model inis for all 16 models (incl. the 6 without test inis) are listed in your brief.

KNOWN FACTS: all 73 test models run on A0. EFTflag=0 vs plain CAMB differ ~2e-5 (reionization setting). CAMB v2 changed tensor Cls ~14% vs the pre-rebase tree even in GR. Known non-converged models: 3_RPH_MB_const_1, 5_Kmimic_2, 5_Horava_SolSyst_*. Python getters returning strings (model_name etc.) are fine on A0. A parallel audit workflow is checking the ctypes layout: if cached values look garbled or shifted, report a layout suspicion rather than a physics bug.

FINDINGS: a defect in the code under audit. Severity: S1 = wrong physics in a supported model, or crash; S2 = wrong only in an edge case / parameter region, or silently poor accuracy; S3 = robustness, docs, naming, dead code. DISAGREEMENT PROTOCOL: when your derivation and the code disagree, settle it by a third route (Bianchi identity / constraint propagation, a cross-basis or cross-code run, or a known limit) before filing; give both expressions side by side in the finding; if the paper itself is wrong, say so (then the code may be right). Do not fix the repo; propose a fix. fix_class: "physics-changing" for anything altering outputs.
MATRIX: report every quantitative check as a matrix entry {model, check, value (worst residual or deviation), threshold, pass, note}; use "n/a" with a reason when a check does not apply. Use model keys exactly: 007p1-pure-eft, 007p2-rph, 007p3-fR, 007p4-mc-quint, 007p5-OL-gamma, 007p6-OL-alpha, 007p7-shiftsym, 007p8-3p1, 008p0-horndeski, 008p1-horava, 008p2-ade, 008p3-kmouflage, 008p4-quintessence, 008p5-beyond-horndeski, 008p6-scaling-cubic, 008p7-ext-galileon.`

const MATRIX_ITEM = {
  type: 'object',
  properties: {
    model: { type: 'string' }, check: { type: 'string' }, value: { type: 'string' },
    threshold: { type: 'string' }, pass: { type: 'boolean' }, note: { type: 'string' },
  },
  required: ['model', 'check', 'value', 'pass'],
}

const FINDING_ITEM = {
  type: 'object',
  properties: {
    slug: { type: 'string' }, title: { type: 'string' },
    severity: { type: 'string', enum: ['S1', 'S2', 'S3'] },
    file: { type: 'string' }, line: { type: 'integer' },
    other_locations: { type: 'array', items: { type: 'string' } },
    repro_script: { type: 'string' }, repro_command: { type: 'string' },
    observed: { type: 'string' }, expected: { type: 'string' },
    code_expression: { type: 'string' }, derived_expression: { type: 'string' }, third_route: { type: 'string' },
    root_cause: { type: 'string' }, proposed_fix: { type: 'string' },
    fix_class: { type: 'string', enum: ['safe-local', 'physics-changing', 'unknown'] },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
  },
  required: ['slug', 'title', 'severity', 'file', 'observed', 'expected', 'root_cause', 'proposed_fix', 'fix_class', 'confidence'],
}

const RESULT_SCHEMA = {
  type: 'object',
  properties: {
    id: { type: 'string' },
    verdict: { type: 'string', enum: ['PASS', 'FINDINGS', 'INCOMPLETE'] },
    summary: { type: 'string' },
    derive_first_followed: { type: 'boolean' },
    sources: { type: 'array', items: { type: 'string' } },
    scripts: { type: 'array', items: { type: 'string' } },
    commands: { type: 'array', items: { type: 'string' } },
    matrix: { type: 'array', items: MATRIX_ITEM },
    findings: { type: 'array', items: FINDING_ITEM },
    not_covered: { type: 'array', items: { type: 'string' } },
  },
  required: ['id', 'verdict', 'summary', 'derive_first_followed', 'sources', 'scripts', 'commands', 'matrix', 'findings', 'not_covered'],
}

const PREP_SCHEMA = {
  type: 'object',
  properties: {
    models: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          key: { type: 'string' }, model_file: { type: 'string' },
          selection: { type: 'string' },
          inis: { type: 'array', items: { type: 'string' } },
          representative_ini: { type: 'string' },
          runs_ok: { type: 'boolean' }, notes: { type: 'string' },
        },
        required: ['key', 'model_file', 'selection', 'inis', 'representative_ini', 'runs_ok', 'notes'],
      },
    },
    findings: { type: 'array', items: FINDING_ITEM },
    notes: { type: 'string' },
  },
  required: ['models', 'findings', 'notes'],
}

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    reproduced: { type: 'boolean' }, refuted: { type: 'boolean' },
    verdict: { type: 'string', enum: ['CONFIRMED', 'PLAUSIBLE', 'REFUTED'] },
    corrected_severity: { type: 'string', enum: ['S1', 'S2', 'S3', 'not-a-bug'] },
    root_cause_agrees: { type: 'boolean' }, reason: { type: 'string' }, evidence: { type: 'string' },
  },
  required: ['reproduced', 'refuted', 'verdict', 'corrected_severity', 'root_cause_agrees', 'reason'],
}

const MODELS = [
  { key: '007p1-pure-eft', file: '07f_designer_models/007p1_Pure_EFT_std.f90', notes: 'Pure EFT: Omega(a) and gamma_i(a) given directly by parametrizations, plus w_DE(a); c and Lambda reconstructed from the background (designer). Check: the background reproduces the input w_DE(a) to 1e-8 (independent integration of the Friedmann equations with that w_DE and the same matter content incl. massive neutrinos); the reconstructed c and Lambda equal the analytic designer expressions you derive from the EFT background equations with Omega(a); cached Omega, gamma_i and their a-derivatives equal the input parametrizations (cover the PEFT test inis: const, linear, power, exp for Omega and gamma1..6, wDE LCDM/wCDM).' },
  { key: '007p2-rph', file: '07f_designer_models/007p2_RPH.f90', notes: 'RPH (reparametrised Horndeski): four functions of time plus w_DE mapped to EFT functions (check the source: EFTCAMB RPH implementation, e.g. arXiv:1405.3590 section on RPH, Gleyzes et al. 1304.4840). Derive the RPH -> EFT map and its inverse independently; compare cached Omega, gamma_i, c, Lambda (and alphas) for the 3_RPH_* test inis; inverse map recovers the inputs.' },
  { key: '007p3-fR', file: '07f_designer_models/007p3_Designer_fR.f90', notes: 'Designer f(R): derive the designer ODE (Song, Hu & Sawicki astro-ph/0610532; Pogosian & Silvestri 0709.0296; EFTCAMB 1312.5742) and solve it independently with scipy at high accuracy for the same w_DE and B0; compare f_R(a) = Omega(a) and the Compton parameter B(a) to 1e-6, and the EFT map (gamma_i = 0, c and Lambda from f). B0 -> 0 must give LCDM. Derive the standard f(R) quasi-static mu(k,a), Sigma(k,a) and compare with compute_QS_mu/Sigma (EFT_mu1/EFT_mu2 cache fields or the QS outputs) and with the full solution at k >> aH.' },
  { key: '007p4-mc-quint', file: '07f_designer_models/007p4_Designer_mc_quintessence.f90', notes: 'Designer minimally coupled quintessence: Omega = gamma_i = 0; derive c and Lambda from rho_DE, P_DE (EFTCAMB normalisation) and check; the reconstructed scalar field quantities (phi_scf, V_scf, dphi... if filled) reproduce w_DE(a); the perturbations reduce to a c_s^2 = 1 scalar. Test inis 4_Des_mc5e_*.' },
  { key: '007p5-OL-gamma', file: '07f_designer_models/007p5_Omega_Lambda_gamma.f90', notes: 'Omega-Lambda parametrisation of Horndeski with gamma functions: find which functions are inputs (Omega, Lambda?, gamma_i), derive c and the background from the EFT background equations; check round trip of inputs and background consistency (no test ini: use the one in the brief).' },
  { key: '007p6-OL-alpha', file: '07f_designer_models/007p6_Omega_Lambda_alpha.f90', notes: 'Omega-Lambda parametrisation with alpha functions: find the inputs (alpha_B, alpha_K, alpha_M/T, ...); derive the alpha -> EFT (Omega, gamma_i) map independently (1404.3713, Gleyzes et al. 1411.3712, review 1907.03150); check cached alphas equal the inputs (round trip) and cached Omega/gamma_i equal the map.' },
  { key: '007p7-shiftsym', file: '07f_designer_models/007p7_ShiftSym_alphaB.f90', notes: 'Shift-symmetric parametrisation of alpha_B: identify the source (file comments/authors; search arXiv for shift-symmetric Horndeski alpha_B parametrisation in EFTCAMB, e.g. Traykova, Bellini, Ferreira et al.), derive the constraint(s) shift symmetry imposes (e.g. relation between alpha_B and the background / other alphas) and verify they hold on the cached functions; check the EFT map.' },
  { key: '007p8-3p1', file: '07f_designer_models/007p8_3p1_wDE_designer.f90', notes: 'EFT of a full 1+3 theory with arbitrary w_DE (designer): identify the 3+1 formulation and source (file comments; Kase & Tsujikawa 3+1 / Gleyzes et al.); check the background reproduces w_DE; derive the map from the 3+1 Lagrangian functions to EFT functions and compare.' },
  { key: '008p0-horndeski', file: '08f_full_models/008p0_Horndeski.f90', notes: 'Full Horndeski G2..G5(phi, X): sympy derivation of the background equations (Kobayashi, Yamaguchi, Yokoyama 1105.5723; De Felice & Tsujikawa) and of the map to EFT functions and alphas (Gleyzes et al. 1304.4840, 1411.3712; 1404.3713). Check the scalar-field background equation residual on the code solution; compare cached Omega, c, Lambda, gamma_i, alphas for 5_hdsk_TG_* inis. The file fortran/eftcamb/08f_full_models/hdsk_coefficients*.c is generated code: regenerate its expressions from your sympy derivation and compare numerically at random inputs. Papers cited in the file: arXiv:2112.06892, 2403.13096.' },
  { key: '008p1-horava', file: '08f_full_models/008p1_Horava.f90', notes: 'Low-energy Horava gravity (Frusciante, Raveri, Vernieri, Hu, Silvestri arXiv:1508.01787): derive the map from (lambda, xi, eta) to EFT functions; derive the background (LCDM-like with a rescaled cosmological gravitational constant) and compare with cached H(a); derive and check the relation imposed by the solar-system variant (HoravaSolarSystem) and the tensor speed. Test inis 5_Horava_eta_*, 5_Horava_SolSyst_*.' },
  { key: '008p2-ade', file: '08f_full_models/008p2_ADE.f90', notes: 'Acoustic Dark Energy (Lin, Hu, Raveri arXiv:1905.12618): k-essence with constant rest-frame sound speed cs2 and w(a) set by (p, a_c, f_ac, wf). Derive rho(a), w(a) and the EFT map from the paper; check the rest-frame c_s^2 computed from the cached EFT functions equals cs2 at all times; background vs the paper. Test ini 5_ADE_1.' },
  { key: '008p3-kmouflage', file: '08f_full_models/008p3_Kmouflage.f90', notes: 'K-mouflage (arXiv:1509.00611, 1809.09958 cited in the file): solve the background ODE independently; derive the EFT map and compare (5_Kmouflage_*). Determine what the Kmimic inis (5_Kmimic_*) are (K-mouflage mimicking a given background?) and check them too.' },
  { key: '008p4-quintessence', file: '08f_full_models/008p4_quintessence.f90', notes: 'Full-mapping quintessence with potentials (5_quint_axion, cubic, exponential, inverse_power, quadratic): Klein-Gordon residual (conformal time) on the code background solution (phi_scf, dphi, ddphi, V_scf cache fields); EFT map Omega = gamma_i = 0 with c and Lambda from phi\' and V (derive the EFTCAMB normalisation); w(a) and rho_DE vs an independent KG + Friedmann integration with the same initial-condition convention (find how ICs/potential amplitude are set).' },
  { key: '008p5-beyond-horndeski', file: '08f_full_models/008p5_Beyond_Horndeski.f90', notes: 'Beyond Horndeski consistent with GW170817 (GLPV with c_T = 1): alpha_T must vanish to machine precision; identify the Lagrangian from the source (file comments, EFTCAMB papers on beyond Horndeski), derive alpha_H and the EFT functions incl. the beyond-Horndeski operator independently and compare (no test ini: use the one in the brief).' },
  { key: '008p6-scaling-cubic', file: '08f_full_models/008p6_Scaling_Cubic.f90', notes: 'Scaling cubic galileon (5_scg.ini). The code says it follows "the mathematica file", which is not in the repo. Identify the source paper (search arXiv: scaling cubic galileon EFTCAMB, Frusciante, Kase, Tsujikawa et al.), re-derive the background (scaling solution) and the EFT map, compare.' },
  { key: '008p7-ext-galileon', file: '08f_full_models/008p7_Extended_Galileon.f90', notes: 'Extended Galileon (De Felice & Tsujikawa tracker solution; EFTCAMB implementation, see file comments, e.g. Peirone, Frusciante, Renzi et al.): derive the tracker background H(a) analytically (parameter xi from the tracker) and compare with cached H; derive the EFT map and compare (no test ini: use the one in the brief).' },
]

const SHARED = [
  { id: '7a-background-coeffs', prompt: 'For every one of the 16 models (representative ini; more inis per model if time allows), on a dense a-grid after the background/perturbation turn-on: (1) background: rho_DE and P_DE built from the cached Omega, c, Lambda (and derivatives) satisfy the DE continuity equation and both Friedmann equations with the cached adotoa, Hdot, Hdotdot and the cached matter densities/pressures; residual relative to the largest term, threshold 1e-6. (2) coefficients: evaluate the reference expressions for every cached factor (EFTeomF, EFTeomN, EFTeomNdot, EFTeomX, EFTeomXdot, EFTeomY, EFTeomG, EFTeomU, EFTeomL, EFTeomM, EFTeomV, EFTeomVdot, EFTeomQ, EFTpiA1, EFTpiA2, EFTpiB1, EFTpiB2, EFTpiC, EFTpiD1, EFTpiD2, EFTpiE, EFTAT, EFTBT, EFTDT, EFT_kinetic, EFT_gradient, alphaB/K/M/T and dots, Meff2, and the s00/s0i/sii/sij-type coefficient fields if the reference defines them) from the cached Omega, c, Lambda, gamma_i and their derivatives, and compare with the cached values: threshold 1e-8 relative (pure algebra); for k-dependent factors use get_scale_evolution at several k. Report the worst residual per model per factor group in the matrix and every mismatch as a finding (identify the exact term).' },
  { id: '7a-evolved-solution', prompt: 'Using get_scale_evolution on A0 for every one of the 16 models at k in {1e-4, 1e-3, 1e-2, 0.1, 1, 10}/Mpc and a dense a-grid after the perturbation turn-on: plug the evolved variables (etak, z, dz, sigma, sigmadot, pi, pidot, pidotdot or delta_phi..., dgrho, dgq, clxc, clxb, clxg, clxr, vb, dgpnu, h-related fields) into the redundant Einstein equations (those CAMB does not use as constraints) and into the pi equation from the reference; residual relative to the largest term must be < 1e-4 (integration tolerance) at all k and times; find where it fails and why. Also check: (a) the 00 and 0i constraints the code imposes hold exactly; (b) switching EFTCAMB_evolve_delta_phi and EFTCAMB_evolve_metric_h on/off (for models that support them) gives spectra equal to < 1e-4; (c) the evolved pi at early times matches the initial-condition prescription (09_EFTCAMB_IC.f90) and the attractor solution derived from the reference pi equation.' },
  { id: '7a-observables', prompt: 'For every one of the 16 models: (1) the Newtonian-gauge potentials Psi and Phi rebuilt from synchronous-gauge variables plus pi by the gauge transformation (derive it, with EFTCAMB pi convention) equal the cached Psi, Phi, PsiDot, PhiDot; (2) the lensing and ISW source terms (EFTLensing, EFTISW cache fields, and the CMB source functions in equations.f90 output) use the metric Phi+Psi (and its time derivative) with no model-dependent factor; (3) at k >> aH (k = 0.1-1/Mpc) and late times (z < 3), mu and Sigma measured from the full solution (mu = -k^2 Psi / (4 pi G a^2 rho_m Delta_m) with the right normalisation; Sigma from Phi+Psi) match compute_QS_mu/Sigma (EFT_mu1/EFT_mu2 cache fields, check meaning) and your analytic QSA expressions (arXiv:1606.05339) to a few per cent where the QSA holds (justify where it does not: e.g. near the sound horizon of the scalar); (4) the sub-horizon CDM growth from P(k,z) at several z < 5 matches an independent growth ODE with your QSA mu to 1% (CDM-only approximation; quantify baryon/neutrino corrections).' },
]

const CROSS = [
  { id: '7c-fluid', prompt: 'Cross-check against CAMB\'s own dark-energy fluid. Run CAMB with EFTflag=0 and dark_energy_model fluid with a tabulated w(a) (DarkEnergyFluid.set_w_a_table, cs2=1 rest-frame) on A0 (and B0) and compare Cls (ell<=2500, TE normalised) and P(k) with: Pure EFT with only w_DE (Omega = gamma_i = 0; e.g. 2_PEFT_wDE_wCDM_1 and a CPL variant) using the same w(a); designer mc quintessence (4_Des_mc5e_*) using its w(a); full-mapping quintessence (5_quint_*) using the w(a) extracted from its cached background (rho_DE, P_DE); ADE (5_ADE_1) using its w(a) and cs2. Use only w > -1 (no phantom crossing). Target 1e-4; make sure the cosmology (H0, densities, neutrinos, DE density today) matches exactly and the w table is dense enough (check by refining). Any deviation above target: localise (background vs perturbations: compare H(z), then transfer functions / get_scale_evolution of clxc) and report.' },
  { id: '7c-bases', prompt: 'Same theory in different bases. Take a full model\'s EFT functions (Horndeski 5_hdsk_TG_1, Extended Galileon (ini in brief), K-mouflage 5_Kmouflage_1) on a dense a-grid (get_eft_functions) and feed them as spline parametrizations (fortran/eftcamb/04f_parametrizations_1D/04p016_spline_parametrizations_1D.f90 or 04p019_spline5_parametrizations_1D.f90: find how nodes are given in the ini) into (a) the gamma basis (Pure EFT 007p1, and Omega_Lambda_gamma 007p5) and (b) the alpha basis (Omega_Lambda_alpha 007p6), with the same w_DE or background. All runs must agree with the original full-model run to ~1e-4 in Cls/P(k) (confirm the residual shrinks when the spline is refined, so it is interpolation error). A disagreement that does not shrink means one basis translation is wrong: localise it by comparing cached EFT functions and coefficient fields between the runs.' },
  { id: '7c-fR-pureEFT', prompt: 'Designer f(R) vs Pure EFT: run designer f(R) (4_Des_F_R_*), extract Omega(a) = f_R(a) on a dense grid and w_DE, then run Pure EFT with a spline Omega(a) = f_R(a), gamma_i = 0 and the same w_DE (find how spline nodes are passed in the ini). Spectra must agree to ~1e-4 (refine the spline to show convergence). Also compare the cached c, Lambda and the pi-equation coefficients between the two runs.' },
]

function verifyAll(res, owner) {
  if (!res) return null
  const fs = res.findings || []
  log(`${owner}: ${res.verdict || ''} ${fs.length} findings`)
  return parallel(fs.map(f => () => {
    const lenses = f.severity === 'S3' ? ['reproduce'] : ['reproduce', 'third-route']
    return parallel(lenses.map(lens => () => agent(`${CTX}

You are an ADVERSARIAL VERIFIER of one physics finding from agent ${owner}. Try to refute it; confirm only what you reproduce and cannot explain away. Do not fix anything. Scripts under ${PHYS}/verify/, outputs under ${AUDIT_ROOT}/work/verify-phys/.
${lens === 'reproduce'
  ? 'LENS: reproduce. Run the repro script, then check the claim with your own independent script. Is the code really doing what the finding says (read the cited code)? Is the effect real and of the stated size?'
  : 'LENS: third route. Independently re-derive the disputed expression from the source paper/action (do not trust either the finding or the code), and test it by a third route: Bianchi identity / constraint propagation, a known limit, or a cross-basis / cross-code run. Decide whether the code, the finding\'s derivation, or the paper is wrong, and whether severity is right.'}

FINDING:
${JSON.stringify(f, null, 2)}

Return reproduced, refuted, verdict (CONFIRMED / PLAUSIBLE / REFUTED), corrected_severity, root_cause_agrees, reason, evidence (numbers).`, {
      label: `verify:${owner}:${f.slug}:${lens}`, phase: 'Verify', schema: VERDICT_SCHEMA,
      effort: f.severity === 'S3' ? 'medium' : undefined,
    }))).then(vs => ({ ...f, owner, verdicts: vs.filter(Boolean) }))
  })).then(v => ({ ...res, findings: v }))
}

phase('Prep')
const prep = await agent(`${CTX}

PREP TASK (id: prep). Produce a validated ini set for all 16 models so other agents can start immediately.
For each model key (list in MATRIX rules) give: the model file under fortran/eftcamb/, the selection (EFTflag and the model selection flag value, read from fortran/eftcamb/09_EFTCAMB_main.f90 allocation logic and the model's read_model_selection), the existing test inis that use it (fortran/eftcamb_test/parameters/*.ini; grep EFTflag / model flags), and one representative ini. For the 6 models without a test ini (007p5 Omega_Lambda_gamma, 007p6 Omega_Lambda_alpha, 007p7 ShiftSym_alphaB, 007p8 3p1 wDE designer, 008p5 Beyond Horndeski, 008p7 Extended Galileon) write inis under fortran/eftcamb_test/audit/inis/<key>.ini (reuse files another agent already put there if they work): DEFAULT(/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/base_params.ini) (or hdsk_base_params.ini when appropriate), keys from the model's read_model_selection / init_model_parameters_from_file, physically sensible parameters (from the source paper if possible) that run on A0 and pass the stability checks. Also write a second ini per new model with clearly non-GR parameters if the first is close to GR. Verify each with harness.run_model on A0 (threads=2) and report runs_ok. Report any defect you hit (e.g. a model that cannot be configured from an ini, wrong key names, crash) as a finding. Write the mapping also to ${PHYS}/models.json.`, { label: 'prep:inis', phase: 'Prep', schema: PREP_SCHEMA })

const modelTable = JSON.stringify(prep ? prep.models : [], null, 1)
const BRIEF = `MODEL INIS (from the prep agent; paths relative to the repo unless absolute):\n${modelTable}`

const referenceBranch = async () => {
  const refPrompt = (which) => `${CTX}

${BRIEF}

REFERENCE DERIVATION ${which} (id: 7a-ref-${which}). You are one of two agents deriving the shared EFT equations completely independently (you will not see the other derivation). Do NOT read fortran/eftcamb/06_abstract_EFTCAMB_model.f90 or any EFT model file before your module is written.
Derive, in sympy, from the EFTCAMB unitary-gauge EFT action (conventions of arXiv:1312.5742 and numerical notes arXiv:1405.3590: Omega, Lambda, c, and the second-order operators with EFTCAMB's dimensionless gamma_1..gamma_6 normalisations, as in the papers):
 (1) background equations (both Friedmann equations; c and Lambda in terms of Omega, H, matter);
 (2) the linear scalar perturbation equations in synchronous gauge in CAMB variables (conformal time, adotoa, k, etak, z = h'/(2k)?, sigma, dgrho, dgq, dgpi... take the GR conventions from CAMB's fortran/equations.f90 GR code and the CAMB notes): all four Einstein equations and the pi (Stuckelberg) equation;
 (3) the tensor equation;
 (4) the no-ghost and gradient-stability conditions (scalar and tensor);
 (5) the alpha-basis functions (Bellini-Sawicki alphas and M*^2) in terms of Omega and gamma_i;
 (6) the quasi-static mu(k,a), Sigma(k,a) (arXiv:1606.05339 conventions).
Then map your equations onto the cached factor names defined in the numerical notes (EFTeomF, EFTeomN, EFTeomNdot, EFTeomX, EFTeomXdot, EFTeomY, EFTeomG, EFTeomU, EFTeomL, EFTeomM, EFTeomV, EFTeomVdot, EFTeomQ, EFTpiA1, EFTpiA2, EFTpiB1, EFTpiB2, EFTpiC, EFTpiD1, EFTpiD2, EFTpiE, EFTAT, EFTBT, EFTDT, EFT_kinetic, EFT_gradient, EFT_mu1, EFT_mu2, s00..sij fields if defined): take the DEFINITIONS from the paper (latest arXiv version), and for factors not defined there you may read where equations.f90 USES them to learn which equation term they multiply, then derive them yourself.
 (7) Show in sympy that the Bianchi identity holds: the conformal-time derivative of the 00 and 0i constraints vanishes on the remaining Einstein equations + pi equation + matter conservation (symbolic, with background functions as undefined functions and background equations substituted). If it fails, find whether your equations or the paper's are wrong.
Deliverable: ${PHYS}/7a-ref-${which}/eft_reference_${which}.py exposing plain numpy functions taking a dict of EFTCAMB_timestep_cache field values (cache field names; arrays allowed) and k: background_residuals(c), factors(c, k) -> dict of predicted cached factors, einstein_residuals(c, k), pi_residual(c, k), tensor_coeffs(c), alphas(c), stability(c, k), qsa(c, k); plus ${PHYS}/7a-ref-${which}/bianchi.py. Document every equation with its source. Sanity test your module on one A0 run (e.g. designer f(R) 4_Des_F_R_wCDM_1 and a Horndeski model) and record the residuals in the matrix (check names prefixed "ref-${which}:"), but do NOT treat disagreements with the code as findings yet unless a third route confirms them.`
  const [ra, rb] = await parallel([
    () => agent(refPrompt('A'), { label: '7a-ref-A', phase: 'Reference', schema: RESULT_SCHEMA }),
    () => agent(refPrompt('B'), { label: '7a-ref-B', phase: 'Reference', schema: RESULT_SCHEMA }),
  ])
  const rec = await agent(`${CTX}

${BRIEF}

RECONCILIATION (id: 7a-ref-final). Two agents derived the shared EFT equations independently:
 A: ${PHYS}/7a-ref-A/ (summary: ${JSON.stringify(ra && ra.summary)})
 B: ${PHYS}/7a-ref-B/ (summary: ${JSON.stringify(rb && rb.summary)})
Compare A and B numerically for every function (background residuals, every predicted factor, residual functions, tensor, alphas, stability, qsa) on random inputs and on cache values from several A0 runs. For every disagreement, decide which is right by re-deriving that term yourself from the papers and by the Bianchi identity (run both bianchi.py; fix or redo if one fails). Do not consult the EFT Fortran expressions to decide. Produce the final module ${PHYS}/eft_reference.py (same API: background_residuals, factors, einstein_residuals, pi_residual, tensor_coeffs, alphas, stability, qsa) and ${PHYS}/bianchi.py, with docstrings citing equations, and a list of every factor and whether A and B agreed. Then (now allowed) compare the final reference factors with the code's cached values for 4 models (designer f(R), a PEFT gamma model, Horndeski, K-mouflage) and list mismatching factors (as findings only if the Bianchi check/third route supports the reference). Matrix: one entry per model per factor group.`, { label: '7a-ref-final', phase: 'Reference', schema: RESULT_SCHEMA })
  const shared = await parallel(SHARED.map(s => () => agent(`${CTX}

${BRIEF}

SHARED-EQUATION CHECK (id: ${s.id}). Use the reconciled reference module ${PHYS}/eft_reference.py (import it; read its docstrings; reconciliation summary: ${JSON.stringify(rec && rec.summary)}). Do not trust the code's equations; trust the reference only where its Bianchi check passed (see ${PHYS}/bianchi.py).
${s.prompt}`, { label: s.id, phase: 'Shared', schema: RESULT_SCHEMA }).then(r => verifyAll(r, s.id))))
  return { ref_A: ra, ref_B: rb, ref_final: await verifyAll(rec, '7a-ref-final'), shared }
}

const modelBranch = () => parallel(MODELS.map(m => () => agent(`${CTX}

${BRIEF}

PER-MODEL DERIVE-FIRST CHECK (id: ${m.key}; code: fortran/eftcamb/${m.file}).
Step 1 (before opening the model's Fortran, its read/feedback routines excepted for parameter names): identify the theory and the source paper(s); write a sympy derivation of the model -> EFT-functions map (Omega, c, Lambda, gamma_i and/or alphas, plus w_DE / background) and the model-specific background equations, to ${PHYS}/${m.key}/derivation.py. Step 2: run the model on A0 (representative ini and the other inis of this model), evaluate your expressions on the code's background/parameters and compare with the cached EFT functions (get_eft_functions on a dense grid), threshold 1e-6 unless stated. Step 3: model-specific checks below. Step 4: only now read the model's Fortran to localise any disagreement, and apply the disagreement protocol.
Model notes: ${m.notes}
Matrix check names: "map:<function>", "background:<check>", "specific:<check>".`, { label: `7b:${m.key}`, phase: 'Models', schema: RESULT_SCHEMA }).then(r => verifyAll(r, m.key))))

const crossBranch = () => parallel(CROSS.map(c => () => agent(`${CTX}

${BRIEF}

CROSS-CHECK (id: ${c.id}). ${c.prompt}
Matrix: one entry per model involved ("cross:<what>").`, { label: c.id, phase: 'Cross', schema: RESULT_SCHEMA }).then(r => verifyAll(r, c.id))))

const [reference, models, cross] = await parallel([referenceBranch, modelBranch, crossBranch])
return { prep, reference, models, cross }
