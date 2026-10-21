# EFTCAMB audit: final regression (HEAD vs pre-fix baseline)

HEAD = `eftcamb-audit` at 4819c8f (all approved fixes, A9 included); updated for the Q3/Q11/Q12 commits (eff41e9), see the
last section. The later validation fix 5b214ea (turn-on
before tau = 1e-12 Mpc refused instead of crashing) changes no shipped model (8 models rechecked bit-identical). Baseline A0 = the pre-fix code (f6f4ba9) as frozen
in the audit tree `trees/A0`: `runs/A0` for the 73 shipped test models, plus a fresh A0 run of the 12 audit inis
(`inis/0*.ini`, without the `*.ini.ini` byproducts) on the same tree (`runs/A0aud`). F2 = the integration state
before A9 (`runs/F2`). All runs through `checks/harness.py` (`run_model`, 1 thread, 4 jobs), default accuracy.
The inis themselves changed between A0 and HEAD only by C4 (c) (`EFT_mass_stability = F` in 14 inis) and C8
(`EFTCAMB_stability_threshold = 1.e-8`); every run uses the ini of its own checkout.

Columns: max fractional difference HEAD vs A0 for l = 2..2500 of the lensed total C_l (TE normalised by
sqrt(TT EE)), C_l^phiphi (PP), linear P(k) (all z, all k), sigma8, H(z) on 0 <= z <= 1e4. Runtime = ini read +
`get_results` in the worker (1 thread, 4 concurrent jobs; noisy at the 10% level). "vs F2" = HEAD vs F2 (worst C_l
or P(k); "bit" = bit-identical). The last column gives the worst change (max of C_l, P(k), sigma8) of the same model in
the per-agent tables measured against A0 while the fixes were made: core = A/B/C fixes in equations/results/camb
(`fixA/runs/groupC9/table_vs_A0.txt`), stab = stability/IC/RGR group (`agentB/final_vs_A0.txt`), models = model files
and parametrizations (`me/gate_table.txt`); "-" = model not in that table.

## Summary

- 85 inis (73 test models + 12 audit inis). HEAD vs F2: every run bit-identical (A9, the stop -> GlobalError
  refactor, changes no result), with the same status.
- HEAD vs A0: 1 bit-identical (1_EFT_GR, EFTflag = 0), 4 within 1e-5 or below the C_l noise (2_PEFT_gamma4_lin_1,
  gamma4_power_1/2, wDE_LCDM_1: 1.5e-10 TT, from D F-136, doubled thermal-history sampling), 76 changed
  above 1e-5, 2 status changes (008p7-ext-galileon(_2): error -> ok, B15/B16), 2 refused in both (007p8-3p1(_2): A0
  "background solver failed" (A4), HEAD "theory unstable", see open questions in REVIEW.md).
- Common floor: most models with active EFT perturbations move by TT ~2e-5, EE ~1.1e-5, P(k) <= 1e-6. This comes
  from C9 (integration stops exactly at the EFT turn-on, extra source samples around it): in the core group the
  change appears entirely at the C9 commit (`fixA/runs/groupC9/table_vs_groupD.txt`; groupD, which already has B3,
  is within 8e-9 of the previous state for these models). A near-GR check (PEFT Omega = 1e-9 or Gamma1 = 1e-9,
  RGR skipped) vs EFTflag = 0: TT 2.7e-6 (A0) -> 2.3e-5 (HEAD) with turn-on 0.01, 2.6e-6 with turn-on 1e-8, i.e.
  the GR-limit consistency is now at the 2e-5 level (CAMB's own time-sampling sensitivity is ~6e-5).
- Runtime: median HEAD/A0 0.98; total 2198 s -> 2268 s (+3%). Larger: 4_Des_F_R_LCDM_2 39.5 -> 109.5 s (C3, earlier
  turn-on 0.093 -> 0.016 where the f(R) pi equation is stiff; after Q3 40.1 s, see the update below), 008p5-beyond-horndeski(_2) 17.9/12.9 -> 28.3/19.1 s
  (B16, background from the model's own solution), 5_scg 20.2 -> 25.2 s (B15, massive nu in the ODE),
  007p5/007p6 +12-18%, 008p7-ext-galileon(_2) 204/177 s (failed after 58/31 s in A0).

## Changed models and responsible items

| model(s) | biggest change HEAD vs A0 | items |
|---|---|---|
| 2_PEFT_Omega_const_1/2, gamma3_*, gamma5_const_1, gamma5_exp_1, gamma4_const_1 | TT 1.2e-2 (Omega_const_1, gamma5_const_1), 5-8e-3 (gamma3), P(k) <= 7e-5, sigma8 <= 2.5e-5 | C9 (turn-on jump integrated exactly and sampled; removes the A0 time-sampling error of 1-5e-2, F-029) |
| other 2_PEFT_* with EFT on, 4_Des_mc5e_*, 4_Des_F_R_LCDM_1, wCDM_1/2, 007p7-shiftsym(_2) | TT ~2e-5, EE 1.1e-5 (C9 floor); f(R) P(k) 2-4e-6 | C9; f(R) also stab group (C3 RGR, <= 3e-6) and B6 (Hdotdotdot, get_eft_functions only) |
| 3_RPH_MB_const_1 | P(k) x100 at a single k (the IC pole), TT 1.0, sigma8 5.7e-6 | C2 (pi IC regularisation; the A0 spike and zeros across k_c are gone, default and boost2 agree); C9 TT 7e-2 |
| 3_RPH_MB_OmegaDE_1, MK_const_1, KT/MK_OmegaDE_1, MK_power_1, KT_power_1 | TT 2.7e-3 / 6.5e-3 / 1e-4 | C9; stab group (C2/C8) <= 4.5e-5 |
| 4_Des_F_R_LCDM_2 | P(k) 0.18 (k > 3/Mpc; 1.6e-2 at k <= 1/Mpc), PP 3.6e-3, sigma8 2.3e-4 | C3 + Q3 (dimensionless, observable-based return to GR at k = 1/Mpc: turn-on 0.0931 -> 0.0393; A0 is the less converged one: vs an early turn-on A0 P(k <= 1/Mpc) 1.6e-2, HEAD 4e-4); runtime 1.08x (run alone) |
| 5_ADE_1 | 8.6e-5 C_l, P(k) 7.3e-5, H 3.2e-5 | B21 (10x background tables), C3 (turn-on 1e-10 -> 1.16e-6), C9 |
| 5_Horava_* | TT 0.49 (SolSyst_3), 0.38 (eta_4), 4.3e-2 (SolSyst_2); sigma8 <= 5.7e-3 | B7 (Hubble derivatives from the model's H(a), arXiv:1508.01787), B3 (SolSyst_3 TT 0.40 and eta_4 5e-3 at the B3 step of the core group) and C9 (0.11 / 3e-2) |
| 5_Kmimic_1/2/3 | P(k) 18 / 530 / 0.86 (single-k IC poles removed), TT 0.94 / 0.92 / 0.22, sigma8 0.776 -> 0.616 etc. (2.1e-1) | C2 (IC regularisation, k-independent IC time), B1 (EFTpiDfunction), C9 (core group 0.33 TT, P(k) x9), B21 (10000-point table, c6b0fb0) |
| 5_Kmouflage_1/2/3 | TT 1.3e-2 / 7.7e-2 / 4e-4; H 7-8e-5 | C9; B21 (tables: 2.4e-4 C_l, 7e-5 H) |
| 5_hdsk_TG_1/2/3 | sigma8 0.14 / 0.14 / 0.022; C_l <= 2e-4 | C7 (EFT background H0 adopted); TG_2 also C2 (P(k) 7e-4, now agrees with TG_1 delta_phi) |
| 5_quint_* | C_l 1-2.4e-4, H 2-6e-5; quint_cubic P(k) 4.8e-4 | B21 (tables); C3 (turn-on 1e-8 -> 1-2.5e-3, frozen field, <= 3e-5); quint_cubic P(k) from the stab group |
| 5_scg | TT 0.21, P(k) 1.4, sigma8 7.2e-2, H 1e-2 | B15 + B8 + B23 (massive nu in the background, third derivatives, _dl literal) |
| 007p5-OL-gamma(_2), 007p6-OL-alpha(_2) (audit inis) | TT 1.5e-2 / 4.6e-3 / 4.9e-3 / 3.4e-3, P(k) <= 5.4e-4 | C9 and/or B3 (the tree without both, `fixA/attr/noB3`, is A0 bit-for-bit; not separated further) |
| 008p5-beyond-horndeski(_2) (audit inis) | H(z) 0.86 (A0: CAMB H(z=0) = 0.52 H0), C_l O(1) | B16 (EFT functions and H from the model's own background) |
| 008p7-ext-galileon(_2) (audit inis) | A0 error -> HEAD ok | B15 + B16 (tracker background at all a, massive nu) |

## EFTflag = 0 vs plain CAMB (trees/B0)

1_EFT_GR.ini on HEAD vs plain CAMB v2 (`trees/B0`), max fractional difference over every C_l type, P(k), sigma8,
H(z):

| case | Reion.max_redshift = 50 on both | default (EFTCAMB 1000, CAMB 50) |
|---|---|---|
| flat (ini as shipped, with source windows) | C_l <= 1.4e-13 (BB), PP 2.1e-13, P(k)/sigma8/H identical | EE 2.3e-5, BB 1.6e-5, TT 1.1e-6 |
| flat, no source windows | <= 5.7e-14 | EE 2.3e-5 |
| omk = -0.05, no source windows | <= 1.3e-13 (BB), TT 3.2e-14 | EE 1.9e-5, BB 1.3e-5 |

With max_redshift = 50 on both, the residual is roundoff (<= 1.4e-13), also with curvature (B9 restored Kf(2);
A0 had TT 6.4e-3 at l = 2 for omk = -0.05). The source windows were removed for omk != 0 because CAMB refuses
windows in non-flat models.

## Per-model table

| model | status A0 -> HEAD | TT | EE | TE | BB | PP | P(k) | sigma8 | H(z) | t A0 / HEAD [s] | vs F2 | agent tables (core / stab / models) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1_EFT_GR | ok (bit-identical) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 9.2 / 7.9 | bit | 0 / 0 / 0 |
| 2_PEFT_Omega_const_1 | ok | 1.2e-02 | 1.1e-05 | 7.1e-03 | 5.2e-04 | 7.1e-04 | 1.2e-06 | 4.4e-07 | 0 | 27.0 / 24.5 | bit | 1.3e-02 / 1.1e-06 / 0 |
| 2_PEFT_Omega_const_2 | ok | 3.3e-03 | 1.1e-05 | 2.4e-03 | 1.7e-05 | 2.5e-05 | 4.9e-07 | 8.2e-08 | 0 | 26.8 / 24.4 | bit | 4.3e-03 / 3.9e-07 / 0 |
| 2_PEFT_Omega_exp_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.9e-07 | 5.5e-07 | 4.5e-07 | 1.4e-11 | 0 | 30.3 / 28.1 | bit | 2.3e-05 / 9.7e-08 / 2.7e-13 |
| 2_PEFT_Omega_lin_1 | ok | 4.3e-05 | 1.1e-05 | 4.6e-05 | 1.8e-06 | 3.0e-06 | 2.7e-07 | 4.5e-12 | 0 | 26.8 / 24.2 | bit | 6.2e-05 / 1.3e-09 / 0 |
| 2_PEFT_Omega_lin_2 | ok | 2.0e-05 | 1.1e-05 | 1.2e-05 | 2.3e-07 | 3.2e-07 | 2.3e-07 | 2.2e-11 | 0 | 26.6 / 25.4 | bit | 2.1e-05 / 5.0e-10 / 0 |
| 2_PEFT_Omega_power_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.9e-07 | 5.5e-07 | 2.8e-07 | 3.1e-11 | 0 | 27.8 / 26.8 | bit | 2.3e-05 / 9.9e-08 / 0 |
| 2_PEFT_Omega_power_2 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.5e-07 | 4.6e-07 | 5.4e-07 | 7.7e-11 | 0 | 27.9 / 26.9 | bit | 2.3e-05 / 8.9e-09 / 0 |
| 2_PEFT_gamma1_const_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.4e-07 | 4.5e-07 | 2.8e-07 | 3.9e-11 | 0 | 13.7 / 12.1 | bit | 2.3e-05 / 0 / 0 |
| 2_PEFT_gamma1_const_2 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.4e-07 | 4.5e-07 | 2.8e-07 | 3.9e-11 | 0 | 13.4 / 12.0 | bit | 2.3e-05 / 0 / 0 |
| 2_PEFT_gamma3_const_1 | ok | 7.7e-03 | 1.8e-05 | 3.9e-03 | 2.3e-04 | 3.0e-04 | 6.9e-05 | 2.5e-05 | 0 | 14.5 / 12.4 | bit | 7.8e-03 / 5.3e-13 / 0 |
| 2_PEFT_gamma3_const_2 | ok | 6.8e-03 | 1.3e-05 | 3.6e-03 | 4.6e-05 | 5.5e-05 | 2.1e-05 | 1.0e-05 | 0 | 14.7 / 12.2 | bit | 7.1e-03 / 2.7e-14 / 0 |
| 2_PEFT_gamma3_exp_1 | ok | 6.3e-03 | 1.9e-05 | 3.4e-03 | 2.4e-04 | 3.1e-04 | 5.9e-05 | 2.3e-05 | 0 | 16.6 / 14.7 | bit | 6.4e-03 / 2.7e-13 / 0 |
| 2_PEFT_gamma3_exp_2 | ok | 5.3e-03 | 1.1e-05 | 3.0e-03 | 3.2e-05 | 4.1e-05 | 1.5e-05 | 7.1e-06 | 0 | 16.3 / 14.4 | bit | 5.6e-03 / 3.7e-14 / 0 |
| 2_PEFT_gamma3_lin_1 | ok | 3.1e-04 | 1.1e-05 | 3.0e-04 | 3.2e-06 | 4.1e-06 | 1.6e-06 | 6.6e-07 | 0 | 14.1 / 12.3 | bit | 3.8e-04 / 3.9e-14 / 0 |
| 2_PEFT_gamma3_lin_2 | ok | 1.8e-04 | 1.1e-05 | 1.0e-04 | 2.4e-07 | 3.9e-07 | 5.2e-07 | 1.4e-07 | 0 | 14.1 / 12.1 | bit | 2.4e-04 / 2.6e-14 / 0 |
| 2_PEFT_gamma3_power_1 | ok | 6.4e-03 | 2.1e-05 | 3.5e-03 | 1.8e-04 | 2.2e-04 | 5.2e-05 | 2.1e-05 | 0 | 14.5 / 14.1 | bit | 6.5e-03 / 5.2e-13 / 0 |
| 2_PEFT_gamma3_power_2 | ok | 5.3e-03 | 1.1e-05 | 3.0e-03 | 3.1e-05 | 3.9e-05 | 1.4e-05 | 6.8e-06 | 0 | 13.9 / 13.2 | bit | 5.6e-03 / 2.4e-14 / 0 |
| 2_PEFT_gamma4_const_1 | ok | 8.9e-04 | 1.1e-05 | 1.2e-05 | 1.0e-05 | 2.9e-05 | 4.0e-05 | 2.7e-06 | 0 | 23.0 / 22.6 | bit | 9.3e-04 / 7.1e-14 / 0 |
| 2_PEFT_gamma4_lin_1 | ok | 1.5e-10 | 4.4e-09 | 2.6e-09 | 1.6e-11 | 4.7e-12 | 0 | 0 | 0 | 13.7 / 11.7 | bit | 8.0e-09 / 0 / 0 |
| 2_PEFT_gamma4_power_1 | ok | 1.5e-10 | 4.4e-09 | 2.6e-09 | 1.6e-11 | 4.7e-12 | 0 | 0 | 0 | 14.4 / 12.7 | bit | 8.0e-09 / 0 / 0 |
| 2_PEFT_gamma4_power_2 | ok | 1.5e-10 | 4.4e-09 | 2.6e-09 | 1.6e-11 | 4.7e-12 | 0 | 0 | 0 | 14.5 / 12.8 | bit | 8.0e-09 / 0 / 0 |
| 2_PEFT_gamma5_const_1 | ok | 1.2e-02 | 2.2e-05 | 5.7e-03 | 5.4e-05 | 8.1e-05 | 1.2e-06 | 3.8e-09 | 0 | 14.5 / 13.0 | bit | 1.2e-02 / 7.0e-14 / 0 |
| 2_PEFT_gamma5_exp_1 | ok | 8.7e-03 | 2.0e-05 | 4.7e-03 | 3.7e-05 | 5.4e-05 | 7.0e-07 | 4.0e-09 | 0 | 14.9 / 14.0 | bit | 9.1e-03 / 5.8e-14 / 0 |
| 2_PEFT_gamma5_exp_2 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.8e-07 | 4.9e-07 | 2.0e-07 | 1.0e-11 | 0 | 14.3 / 13.5 | bit | 2.3e-05 / 2.4e-14 / 0 |
| 2_PEFT_gamma6_const_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.4e-07 | 4.5e-07 | 2.8e-07 | 3.9e-11 | 0 | 13.3 / 13.1 | bit | 2.3e-05 / 0 / 0 |
| 2_PEFT_gamma6_exp_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.4e-07 | 4.5e-07 | 2.8e-07 | 3.9e-11 | 0 | 14.3 / 13.1 | bit | 2.3e-05 / 0 / 0 |
| 2_PEFT_gamma6_lin_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.4e-07 | 4.5e-07 | 2.8e-07 | 3.9e-11 | 0 | 13.9 / 13.1 | bit | 2.3e-05 / 0 / 0 |
| 2_PEFT_gamma6_power_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.4e-07 | 4.5e-07 | 2.8e-07 | 3.9e-11 | 0 | 14.1 / 13.2 | bit | 2.3e-05 / 0 / 0 |
| 2_PEFT_wDE_LCDM_1 | ok | 1.5e-10 | 4.4e-09 | 2.6e-09 | 1.6e-11 | 4.7e-12 | 0 | 0 | 0 | 13.6 / 12.8 | bit | 8.0e-09 / 0 / 0 |
| 2_PEFT_wDE_wCDM_1 | ok | 1.5e-05 | 7.7e-06 | 7.9e-06 | 3.3e-07 | 4.4e-07 | 3.2e-07 | 1.7e-10 | 0 | 26.2 / 27.3 | bit | 1.7e-05 / 1.3e-08 / 0 |
| 3_RPH_KT_OmegaDE_1 | ok | 1.9e-05 | 5.6e-05 | 2.7e-05 | 1.7e-07 | 3.8e-06 | 1.5e-06 | 9.5e-10 | 0 | 37.9 / 38.4 | bit | 1.1e-04 / 2.6e-06 / 0 |
| 3_RPH_KT_power_1 | ok | 2.2e-05 | 1.2e-05 | 1.2e-05 | 3.2e-07 | 4.5e-07 | 3.1e-07 | 6.0e-10 | 0 | 14.1 / 13.7 | bit | 2.3e-05 / 1.4e-09 / 0 |
| 3_RPH_MB_OmegaDE_1 | ok | 2.7e-03 | 5.3e-05 | 2.1e-04 | 3.5e-05 | 3.7e-04 | 1.7e-05 | 4.4e-09 | 0 | 80.1 / 82.5 | bit | 3.3e-03 / 4.5e-05 / 0 |
| 3_RPH_MB_const_1 | ok | 1.0e+00 | 1.0e+00 | 9.5e-01 | 2.1e-02 | 2.3e-01 | 1.0e+02 | 5.7e-06 | 0 | 37.7 / 38.2 | bit | 7.0e-02 / 1.0e+02 / 0 |
| 3_RPH_MK_OmegaDE_1 | ok | 9.0e-05 | 5.6e-05 | 2.0e-05 | 1.4e-06 | 1.8e-05 | 2.7e-06 | 1.1e-09 | 0 | 59.0 / 60.4 | bit | 1.1e-04 / 1.3e-06 / 0 |
| 3_RPH_MK_const_1 | ok | 6.5e-03 | 1.2e-05 | 3.9e-03 | 1.7e-05 | 2.4e-05 | 4.9e-07 | 2.3e-08 | 0 | 36.3 / 38.6 | bit | 8.1e-03 / 3.3e-07 / 0 |
| 3_RPH_MK_power_1 | ok | 7.1e-05 | 1.2e-05 | 4.5e-05 | 4.4e-07 | 4.7e-05 | 3.4e-07 | 1.7e-09 | 0 | 44.9 / 45.9 | bit | 8.7e-05 / 6.5e-07 / 0 |
| 4_Des_F_R_LCDM_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 2.2e-07 | 2.8e-07 | 4.3e-06 | 4.1e-09 | 0 | 27.1 / 25.4 | bit | 2.3e-05 / 3.3e-06 / 0 |
| 4_Des_F_R_LCDM_2 | ok | 3.0e-05 | 8.4e-06 | 1.3e-05 | 1.1e-03 | 3.6e-03 | 1.8e-01 | 2.3e-04 | 0 | 39.5 / 40.1 | Q3: P(k) 1.3e-1, PP 1.0e-4, sigma8 5.3e-6 | 1.9e-04 / 1.3e-01 / 0 |
| 4_Des_F_R_wCDM_1 | ok | 2.0e-05 | 1.0e-05 | 1.1e-05 | 2.0e-07 | 2.4e-07 | 3.4e-06 | 8.2e-12 | 0 | 24.2 / 24.1 | bit | 2.2e-05 / 2.4e-06 / 0 |
| 4_Des_F_R_wCDM_2 | ok | 2.0e-05 | 1.0e-05 | 1.1e-05 | 2.0e-07 | 2.3e-07 | 2.6e-06 | 1.2e-09 | 0 | 24.1 / 24.3 | bit | 2.2e-05 / 2.6e-06 / 0 |
| 4_Des_mc5e_CPL_1 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 3.5e-07 | 4.6e-07 | 2.5e-07 | 1.3e-10 | 0 | 29.5 / 31.0 | bit | 2.3e-05 / 1.0e-09 / 0 |
| 4_Des_mc5e_CPL_2 | ok | 1.5e-05 | 7.8e-06 | 8.1e-06 | 3.3e-07 | 4.4e-07 | 3.8e-07 | 1.7e-10 | 0 | 30.1 / 31.4 | bit | 1.7e-05 / 5.0e-09 / 0 |
| 4_Des_mc5e_taylor_1 | ok | 2.0e-05 | 1.0e-05 | 1.1e-05 | 3.6e-07 | 4.7e-07 | 6.2e-07 | 3.3e-10 | 0 | 30.2 / 31.8 | bit | 2.2e-05 / 8.0e-09 / 0 |
| 4_Des_mc5e_taylor_2 | ok | 2.0e-05 | 1.0e-05 | 1.1e-05 | 3.5e-07 | 4.6e-07 | 4.5e-07 | 2.1e-10 | 0 | 30.1 / 31.3 | bit | 2.2e-05 / 1.3e-08 / 0 |
| 4_Des_mc5e_turning_point_1 | ok | 1.5e-05 | 7.7e-06 | 8.0e-06 | 3.3e-07 | 4.4e-07 | 2.0e-07 | 1.0e-10 | 0 | 30.0 / 31.3 | bit | 1.7e-05 / 9.5e-09 / 0 |
| 4_Des_mc5e_turning_point_2 | ok | 1.9e-05 | 9.7e-06 | 1.1e-05 | 3.7e-07 | 5.0e-07 | 3.6e-07 | 1.9e-10 | 0 | 30.7 / 31.5 | bit | 2.1e-05 / 4.9e-08 / 0 |
| 4_Des_mc5e_wCDM_1 | ok | 1.5e-05 | 7.7e-06 | 7.9e-06 | 3.3e-07 | 4.4e-07 | 3.2e-07 | 1.7e-10 | 0 | 29.9 / 30.8 | bit | 1.7e-05 / 1.3e-08 / 0 |
| 4_Des_mc5e_wCDM_2 | ok | 2.1e-05 | 1.0e-05 | 1.1e-05 | 3.5e-07 | 4.6e-07 | 6.1e-07 | 2.6e-10 | 0 | 30.1 / 30.6 | bit | 2.2e-05 / 5.1e-09 / 0 |
| 5_ADE_1 | ok | 4.6e-05 | 7.6e-05 | 5.9e-05 | 4.6e-05 | 4.7e-05 | 7.3e-05 | 2.3e-05 | 3.2e-05 | 28.7 / 25.0 | bit | 3.7e-05 / 3.4e-05 / 8.6e-05 |
| 5_Horava_SolSyst_1 | ok | 4.8e-04 | 1.1e-05 | 4.1e-04 | 1.2e-05 | 6.0e-05 | 3.0e-07 | 6.7e-10 | 0 | 20.4 / 20.5 | bit | 5.2e-04 / 1.2e-07 / 5.9e-04 |
| 5_Horava_SolSyst_2 | ok | 4.3e-02 | 1.1e-05 | 4.1e-03 | 5.9e-05 | 3.5e-04 | 9.7e-07 | 3.3e-07 | 0 | 20.7 / 20.6 | bit | 4.7e-02 / 3.2e-07 / 5.9e-03 |
| 5_Horava_SolSyst_3 | ok | 4.9e-01 | 6.0e-03 | 3.5e-01 | 3.9e-03 | 5.3e-03 | 6.5e-03 | 3.2e-03 | 0 | 52.4 / 52.8 | bit | 4.6e-01 / 8.8e-06 / 6.2e-01 |
| 5_Horava_eta_1 | ok | 3.3e-04 | 1.1e-05 | 3.7e-04 | 6.0e-06 | 2.5e-05 | 2.4e-07 | 4.2e-09 | 0 | 12.9 / 12.3 | bit | 2.5e-05 / 2.6e-07 / 5.1e-04 |
| 5_Horava_eta_2 | ok | 3.3e-03 | 1.1e-05 | 3.7e-03 | 8.3e-06 | 3.4e-05 | 1.2e-06 | 4.6e-07 | 0 | 12.4 / 12.2 | bit | 1.2e-04 / 3.3e-07 / 5.1e-03 |
| 5_Horava_eta_3 | ok | 3.0e-02 | 1.4e-05 | 3.8e-02 | 6.1e-05 | 7.8e-05 | 9.7e-05 | 4.8e-05 | 0 | 14.0 / 13.8 | bit | 1.0e-03 / 5.5e-07 / 5.4e-02 |
| 5_Horava_eta_4 | ok | 3.8e-01 | 1.7e-03 | 4.0e-01 | 7.0e-03 | 7.3e-03 | 1.1e-02 | 5.7e-03 | 0 | 44.7 / 43.3 | bit | 2.7e-02 / 4.4e-06 / 5.6e-01 |
| 5_Kmimic_1 | ok | 9.4e-01 | 8.6e-01 | 4.0e-01 | 9.1e-01 | 7.4e-01 | 1.8e+01 | 2.1e-01 | 6.6e-05 | 13.5 / 12.2 | bit | 8.9e+00 / 2.1e+01 / 0 |
| 5_Kmimic_2 | ok | 9.2e-01 | 8.1e-01 | 3.8e-01 | 8.7e-01 | 7.0e-01 | 5.3e+02 | 1.7e-01 | 6.5e-05 | 13.1 / 12.2 | bit | 7.1e+00 / 7.5e+02 / 0 |
| 5_Kmimic_3 | ok | 2.2e-01 | 9.8e-02 | 6.4e-02 | 1.4e-01 | 1.5e-01 | 8.6e-01 | 1.6e-02 | 6.9e-05 | 13.1 / 12.4 | bit | 4.4e-01 / 8.6e-01 / 0 |
| 5_Kmouflage_1 | ok | 1.3e-02 | 1.5e-04 | 6.9e-03 | 4.9e-05 | 2.7e-05 | 5.0e-05 | 2.4e-05 | 7.1e-05 | 79.0 / 72.7 | bit | 1.5e-02 / 3.1e-07 / 2.4e-04 |
| 5_Kmouflage_2 | ok | 7.7e-02 | 1.1e-03 | 3.2e-02 | 1.1e-04 | 9.6e-05 | 1.7e-04 | 2.7e-05 | 8.2e-05 | 77.1 / 71.4 | bit | 8.2e-02 / 1.6e-06 / 2.6e-04 |
| 5_Kmouflage_3 | ok | 4.0e-04 | 1.6e-04 | 2.4e-04 | 6.4e-05 | 3.4e-05 | 5.2e-05 | 2.5e-05 | 7.2e-05 | 21.2 / 21.2 | bit | 4.9e-04 / 5.1e-10 / 2.4e-04 |
| 5_hdsk_TG_1 | ok | 1.3e-06 | 3.6e-06 | 2.8e-06 | 6.3e-07 | 7.2e-08 | 1.4e-05 | 1.4e-01 | 0 | 76.6 / 70.7 | bit | 1.4e-01 / 0 / 0 |
| 5_hdsk_TG_2 | ok | 1.6e-04 | 2.0e-04 | 1.0e-04 | 3.6e-06 | 2.5e-05 | 6.4e-04 | 1.4e-01 | 0 | 78.6 / 90.1 | bit | 1.4e-01 / 7.2e-04 / 0 |
| 5_hdsk_TG_3 | ok | 1.7e-06 | 2.4e-06 | 1.3e-06 | 4.9e-07 | 7.9e-08 | 1.7e-05 | 2.2e-02 | 0 | 68.0 / 64.7 | bit | 2.2e-02 / 0 / 0 |
| 5_quint_axion | ok | 4.2e-05 | 6.9e-05 | 4.2e-05 | 1.8e-05 | 1.2e-05 | 2.8e-05 | 3.9e-06 | 6.1e-05 | 25.7 / 26.6 | bit | 2.5e-06 / 2.1e-07 / 9.0e-05 |
| 5_quint_cubic | ok | 4.8e-05 | 7.8e-05 | 4.8e-05 | 2.1e-05 | 9.9e-06 | 4.8e-04 | 5.1e-06 | 5.3e-05 | 30.8 / 30.8 | bit | 3.5e-07 / 4.9e-04 / 9.9e-05 |
| 5_quint_exponential | ok | 2.4e-04 | 1.6e-04 | 1.6e-04 | 3.0e-05 | 1.5e-05 | 2.4e-05 | 1.1e-05 | 2.2e-05 | 23.3 / 24.0 | bit | 2.7e-06 / 3.3e-07 / 1.0e-04 |
| 5_quint_inverse_power | ok | 2.1e-04 | 1.5e-04 | 1.4e-04 | 2.7e-05 | 1.4e-05 | 2.2e-05 | 1.0e-05 | 2.3e-05 | 24.1 / 24.0 | bit | 2.5e-06 / 3.0e-07 / 9.4e-05 |
| 5_quint_quadratic | ok | 1.5e-04 | 1.2e-04 | 1.1e-04 | 5.1e-06 | 1.1e-05 | 9.9e-06 | 3.3e-06 | 4.2e-05 | 23.6 / 25.2 | bit | 2.7e-06 / 2.4e-07 / 1.1e-05 |
| 5_scg | ok | 2.1e-01 | 1.5e-01 | 5.3e-01 | 4.9e-02 | 3.0e-01 | 1.4e+00 | 7.2e-02 | 1.0e-02 | 20.2 / 25.2 | bit | 1.5e-06 / 1.4e-05 / 1.4e+00 |
| 007p5-OL-gamma | ok | 1.5e-02 | 3.8e-04 | 6.5e-03 | 5.5e-07 | 7.2e-05 | 4.0e-04 | 2.7e-07 | 0 | 22.2 / 26.2 | bit | - / - / 0 |
| 007p5-OL-gamma_2 | ok | 4.6e-03 | 6.7e-05 | 1.9e-03 | 1.0e-06 | 1.4e-05 | 2.5e-05 | 5.6e-07 | 0 | 21.0 / 24.4 | bit | - / - / 0 |
| 007p6-OL-alpha | ok | 4.9e-03 | 9.0e-04 | 2.3e-03 | 3.4e-07 | 1.0e-04 | 5.3e-04 | 7.6e-09 | 0 | 37.6 / 42.0 | bit | - / - / 0 |
| 007p6-OL-alpha_2 | ok | 3.4e-03 | 5.9e-04 | 1.1e-03 | 5.0e-07 | 1.4e-04 | 5.4e-04 | 9.4e-09 | 0 | 40.0 / 47.0 | bit | - / - / 0 |
| 007p7-shiftsym | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 4.8e-07 | 3.5e-06 | 3.3e-07 | 3.4e-09 | 0 | 16.1 / 16.9 | bit | - / - / 0 |
| 007p7-shiftsym_2 | ok | 2.1e-05 | 1.1e-05 | 1.2e-05 | 1.5e-06 | 2.8e-05 | 4.2e-07 | 2.0e-09 | 0 | 17.9 / 19.4 | bit | - / - / 0 |
| 007p8-3p1 | error:rc=1 | - | - | - | - | - | - | - | - | 0.4 / 0.3 | error:rc=1/error:rc=1 | - / - / - |
| 007p8-3p1_2 | error:rc=1 | - | - | - | - | - | - | - | - | 0.4 / 0.3 | error:rc=1/error:rc=1 | - / - / - |
| 008p5-beyond-horndeski | ok | 7.9e-01 | 1.6e+00 | 6.6e-01 | 4.6e-01 | 2.8e-01 | 3.4e-01 | 1.5e-01 | 8.6e-01 | 17.9 / 28.3 | bit | - / - / 3.0e+02 |
| 008p5-beyond-horndeski_2 | ok | 2.0e+00 | 7.7e+03 | 3.4e+01 | 4.7e-01 | 6.9e-01 | 9.9e-01 | 1.1e-01 | 8.6e-01 | 12.9 / 19.1 | bit | - / - / 7.8e+03 |
| 008p7-ext-galileon | error:rc=1 -> ok | - | - | - | - | - | - | - | - | 57.6 / 204.1 | bit | - / - / - |
| 008p7-ext-galileon_2 | error:rc=1 -> ok | - | - | - | - | - | - | - | - | 31.2 / 177.4 | bit | - / - / - |


Scripts: scratchpad `final/batch.py` (runs), `final/regress.py` (this table), `final/gr_vs_camb.py`, `final/attr_aud.py`.

## Update 2026-10-21: Q3, Q11, Q12 (da3c2d1, af2504f, eff41e9)

**Regression vs F2** (all 85 inis on eff41e9, `run_model`, 1 thread, 4 jobs; scratchpad `q3/regr.py`, runs in
`q3/regr`): 82 bit-identical, 2 refused in both (007p8-3p1(_2)), 1 changed: 4_Des_F_R_LCDM_2 (Q3: P(k) 0.13 at
k > 3/Mpc, C_l^phiphi 1.0e-4, BB 3.0e-5, TT 1.3e-5, sigma8 5.3e-6). Q11 and Q12 change no result.
`camb.tests.camb_test` + `camb.tests.eftcamb_test` (34 tests): OK; the eftcamb_test reference models are unchanged.

**Q3 trade-off** (4_Des_F_R_LCDM_2; scratch tree with the probe k and threshold read from the environment,
`q3/scan.py`, `q3/errs.py`). Runtime = harness elapsed, 1 thread, one run at a time (A0 32.6-33.4 s in six runs on
the same machine state; pre-Q3 HEAD 94.3 s). Errors = max fractional difference vs a converged reference with turn-on
forced to a = 0.007 (336 s; the turn-on 0.01 run differs from it by P(k <= 10) 7.6e-5, C_l 1.2e-5), linear P(k) at
z = 0, lensed TT/EE/TE (TE / sqrt(TT EE)) and C_l^phiphi for l <= 2500. The C_l column has a ~2e-5 floor (the C9
floor; the turn-on 0.01 reference itself is at 1.2e-5).

| probe k [1/Mpc] | threshold | turn-on a | runtime [s] | P(k <= 1) | P(k <= 10) | sigma8 | lensed C_l | C_l^phiphi |
|---|---|---|---|---|---|---|---|---|
| (A0, Omega-based) | - | 0.0931 | 33.2 | 1.6e-2 | 1.5e-1 | 2.3e-4 | 3.3e-5 | 3.7e-3 |
| 10 (pre-Q3 HEAD) | 1e-3 | 0.0162 | 94.3 | 7.9e-6 | 7.7e-4 | 1.1e-7 | 2.2e-5 | 2.3e-6 |
| 10 | 1e-2 | 0.0274 | 50.3 | 8.2e-5 | 7.7e-3 | 1.1e-6 | 2.4e-5 | 2.3e-5 |
| 10 | 3e-2 | 0.0359 | 39.2 | 2.7e-4 | 2.3e-2 | 3.6e-6 | 2.7e-5 | 7.2e-5 |
| 10 | 5e-2 | 0.0408 | 35.0 | 4.7e-4 | 3.8e-2 | 6.4e-6 | 3.1e-5 | 1.3e-4 |
| 10 | 1e-1 | 0.0500 | 30.4 | 1.1e-3 | 7.4e-2 | 1.6e-5 | 3.3e-5 | 3.0e-4 |
| 3 | 1e-3 | 0.0280 | 49.7 | 8.9e-5 | 8.4e-3 | 1.2e-6 | 2.3e-5 | 2.5e-5 |
| 3 | 3e-3 | 0.0360 | 38.9 | 2.7e-4 | 2.4e-2 | 3.7e-6 | 4.7e-5 | 7.3e-5 |
| 1 | 1e-4 | 0.0273 | 50.5 | 8.1e-5 | 7.6e-3 | 1.1e-6 | 2.2e-5 | 2.2e-5 |
| 1 | 3e-4 | 0.0350 | 39.6 | 2.4e-4 | 2.1e-2 | 3.3e-6 | 2.6e-5 | 6.5e-5 |
| **1 (Q3)** | **5e-4** | **0.0393** | **36.0** | **4.0e-4** | **3.3e-2** | **5.4e-6** | **2.8e-5** | **1.1e-4** |
| 1 | 1e-3 | 0.0462 | 32.2 | 8.1e-4 | 5.8e-2 | 1.1e-5 | 5.9e-5 | 2.1e-4 |
| 1 | 3e-3 | 0.0593 | 27.3 | 2.4e-3 | 1.1e-1 | 3.3e-5 | 9.5e-5 | 6.1e-4 |
| 1 | 1e-2 | 0.0635 | 26.2 | 3.2e-3 | 1.3e-1 | 4.4e-5 | 1.2e-4 | 8.1e-4 |

Chosen: k = 1/Mpc (the linear scales the criterion should protect), threshold 5e-4: the cheapest setting with
P(k <= 1/Mpc) < 1e-3 and all C_l (C_l^phiphi included) at <~ 1e-4; 36.0 s vs A0 33.3 s (1.08x, target <= 1.2x).
Threshold 1e-3 at k = 1 is cheaper (32 s) with lensed C_l 5.9e-5 but C_l^phiphi 2.1e-4. Runtime falls with the
turn-on time (the f(R) pi equation is stiff at early times); 0.0635 is the latest possible turn-on (the c/3H^2
criterion takes over). The price is at small scales: P(k) 3.6e-3 at k <= 3/Mpc, 3.3e-2 at k <= 10/Mpc, 0.14 at
k <= 33/Mpc (A0: 9.2e-2 / 0.15 / 0.15).

**Turn-on and runtime of the models C3 changed** (turn-on from a background-only run for both settings; runtime run
alone, 1 thread; "before" = pre-Q3 setting, 1e-3 at k = 10/Mpc):

| model | turn-on before / after | runtime before / after [s] |
|---|---|---|
| 4_Des_F_R_LCDM_1 | 0.01 / 0.01 (requested time) | 22.0 / 22.5 |
| 4_Des_F_R_LCDM_2 | 0.0162 / 0.0393 | 94.3 / 36.0 (A0 33.3) |
| 4_Des_F_R_wCDM_1 | 0.01 / 0.01 | 21.5 / 22.1 |
| 4_Des_F_R_wCDM_2 | 0.01 / 0.01 | 21.3 / 21.5 |
| 5_ADE_1 | 1.16e-6 / 1.16e-6 | 22.4 / 22.4 |
| 5_quint_axion | 1.11e-3 / 1.11e-3 | 24.2 / 24.2 |
| 5_quint_cubic | unchanged | 28.1 / 28.1 |
| 5_quint_exponential | unchanged | 21.4 / 21.4 |
| 5_quint_inverse_power | unchanged | 21.8 / 21.7 |
| 5_quint_quadratic | unchanged | 22.5 / 22.5 |

All other test models and the audit inis keep their turn-on time (all 73 + 12 checked); their results are
bit-identical to F2.

**Q11**: thetastar / cosmomc_theta with 5_hdsk_TG_1 (Horndeski_shooting = F and T) refused on both paths
("thetastar/cosmomc_theta cannot be used with this model: its H0 (85.5373) is set by the EFT background"); PEFT and
f(R) models still solve theta -> H0 on both paths. **Q12**: the parameter cache of `CAMBparams()`, of a new
`EFTCAMB_parameter_cache()` and after `read_ini` is all zero before the first calculation.

## Reproducers on HEAD

Every script in `repro/` run against HEAD (copies with the A0 tree replaced by HEAD; debug builds of HEAD with the
AD/ADD flags for the three runtime-check reproducers; scratchpad `final/run_repros.py`, logs in `final/rplogs/`).
Result: **fixed** = the expected behaviour is now observed; **error by design** = a clear Python exception
(CAMBError/CAMBValueError) as decided; **partial** = much improved but the reproducer's own target is not fully met;
**n/a** = needs a deleted instrumented tree (the item is checked otherwise, as noted); **upstream** = also in plain
CAMB (B0).

| reproducer | item | result on HEAD | number |
|---|---|---|---|
| 21cm-not-refused-in-python-with-eft | C6 | error by design | CAMBError "EFTCAMB is not compatible with 21 cm sources so far" |
| additional-derivs-wrong-cache | B12 | fixed | EFTcdotdot, Gamma1PP, Gamma2PP vs finite differences 1e-4-7e-4 (were 0 or -2 H cdot) |
| auto-lens-accuracy-enables-nonlinear-for-eft | C6 | by design | automatic non-linear lensing kept; Halofit now gets the model (w0, wa) from results.Params (the script prints the input params, still (-1, 0)) |
| beyond-horndeski-hubble | B16 | fixed | H_BH(z=0)/H0 = 1.00000, DE present (was 0.52) |
| cli-setparams-unallocated-eftcamb | A11 | fixed (run by hand) | -O0 -fcheck=all command-line exe of HEAD passes SetParams and computes 1_EFT_GR |
| cosmomc-theta-exact-rs-use-background | B20 | fixed | theta(T)/theta(F) - 1 = -2.5e-8 (was -1e-3) |
| curvature-silently-wrong-with-eft | C1 | error by design | CAMBError "EFTCAMB does not support spatial curvature" |
| derivst-tensor-uninit-grhog | A20 | fixed (run by hand) | tensor-only 2_PEFT_Omega_const_1 with FP traps on the checked build: clean, BB(2..4) = 0.02559, 0.01482, 0.00877 |
| eft-background-densities-bypass | B18 | fixed | H from get_background_densities / H - 1 = 2.2e-16 |
| eft-getters-gr-segfault | A3 | fixed | all getters return (empty) values with EFTflag = 0 and after a failed init |
| eft-kinetic-cross-term-extra-h0 | B2 | fixed | code / analytic = 1.000; ghost model B rejected |
| eft-par-cache-realloc-abort, -crash | A1 | fixed | get_BAO and get_background(results.Params[.copy()]) rc = 0 with the fresh-run H(z) |
| eft-pickle-repr-loses-model, pickle-deepcopy-eftcamb-segfault | A2 | fixed | copy/pickle/deepcopy/repr round trips keep the model and H(z) |
| eft-turnon-jump-cmb-time-sampling | C9 | partial | 2_PEFT_Omega_const_1 low-l TT: windows on 4.6e-2 -> 1.4e-4; TimeStepBoost = 3 1.9e-2 -> 2.1e-3 (GR 6e-5) |
| eftcamb-import-fails-non-gfortran | A18 | fixed | non-gfortran import branch works |
| evolve-metric-h-inconsistent, metric-h-ic-violates-constraint | C5 | error by design | EFTCAMB_evolve_metric_h = T refused (unsupported parameters) |
| ext-galileon-hdotdot | B15/B16 | fixed | Hdotdot vs exact 5.3e-8 |
| extgal-adotoa-stale-grhonu | B15 | fixed | adotoa vs exact with nu 3e-15, no jump at a = 0.2 |
| fourier-min-a0-zero-order | B13 (F-116) | fixed | 0.2 = expected |
| fr-rgr-threshold-high-k-pk | C3 | fixed | threshold 1e-8 vs 1e-12: P(k=10) 8e-4, sigma8 identical (was percent level) |
| fullmap-bg-eft-stale-background | B16 | n/a | instrumented LAG tree deleted; B16 checked by beyond-horndeski-hubble, ext-galileon-hdotdot |
| fullmap-hdotdotdot-wrong, t10a-fullmap-hdotdotdot-formula | B6 | fixed | Hdotdotdot vs d(Hdotdot)/dtau 7e-4 (FD) / 5.4e-8; f(R) 8e-13 |
| get-eft-functions-gr-segfault | A3 | error by design | CAMBValueError "the results have no EFTCAMB model" |
| get-evolution-pi-highk-tolerance | C9 | partial | default vs IntTolBoost 8: pi 31.9 -> 3.0e-2, Psi 14 -> 1.4e-2, pidot 1101 -> 2.6; tight result sampling-independent (2.6e-9) |
| get-evolution-sort-misalign | A5 | fixed | returned a = requested a (1.3e-9), tau order kept |
| hdsk-h0-not-propagated | C7 | fixed | CAMB H0 = EFT H(z=0) = 85.54, sigma8 error 0; theta does not depend on the input H0 (by design: the EFT H0 is adopted) |
| hdsk-num-eq-implicit-save, hdsk-num-eq-save | A19 | fixed | bit-identical after an evolve_hubble = T run |
| horava-hdot-inconsistent | B7 | fixed | second Friedmann residual 1e-16; stored Hdot vs FD 2e-6 |
| ic1-pi-attractor-pole | C2 | fixed | no pole: delta_c(a=1) 25-58 across k_c (was -7982); default and boost2 agree; the step of P/P_GR at k_c (1.37 -> 0.29) is the regularisation branch (open question) |
| ic1-pidot-eftpid-typo | B1 | n/a | FIX tree deleted; code verified (EFTpiDfunction = D1 + k^2 D2) and K-mimic glitch gone (kmimic-pk-single-mode-glitch) |
| ic1-pidot-transient-noisy-coefficients | C2 | fixed | 5_Kmimic_1 P/P_GR - 1 at k = 1 for turn-on 1e-9/1e-8/1e-7: 0.041/0.041/0.034 (was 1e8); TG_2 vs TG_1 5.7e-9 |
| ic23-pi-radiation-series-at-turn-on | B22 | fixed | IC2 vs IC1 bit-identical; IC 3 runs |
| ic23-stale-cache-gammas | A14 | n/a | instrumented env-var tree; A14 checked by ic3-nan and ic23 reproducers |
| ic3-nan-when-omega-gamma3-zero | A14 | fixed | IC 3 = IC 1, no NaN |
| ini-theta-h0-solved-without-eft | B17 | fixed | thetastar recovered to 6e-8 for 4 EFT models |
| initialize-parameters-bad-value-kills-python, stop-kills-python | A9 | error by design | CAMBValueError / CAMBError with the Fortran message, interpreter alive |
| interp-index-overflow-xfinal | A10 | fixed | checked build passes the stability check (the run then stops at the upstream windows + WantCls=F bound error, F-099) |
| kmimic-pk-single-mode-glitch | C2 | fixed | largest single-k log P deviation -8.9e-3 (GR -9.2e-3; was +1.14) |
| mass-stability-verdict-unreliable | C4 | by design | check experimental and ignored (C4 b); mu values unchanged until the C4 (a) re-derivation |
| mass-stability-wrong-derivatives | B12 | fixed | used = model derivatives |
| minkowski-tachpos-lambdadotdot-zero | B12/B6 | fixed | Minkowski verdict with code and FD EFTLambdadotdot agree (2887/2991 both) |
| nonlinear-halofit-w0wa-lcdm-for-eft | C6 | fixed | P_nl/P_lin EFT vs fluid 5.5e-6 (was 3.3e-2) |
| omp-private-ev-printing-uninit | A23 | fixed | 0 modes with EFTCAMB_printing = T |
| output-pidotdot-zero | B3 | fixed | cached Ndot = FD(N) to FD error (adding the formerly missing term gives O(1)) |
| pade-fourth-derivative, spline-*, spline5-*, steplog-d4-and-integral | B13 | fixed | d4/integrals exact to 1e-15; spline w = -1 identical to EFTwDE = 0; spline order 3 runs |
| param-name-truncation | B14 | fixed | long keys read correctly |
| pert-turn-on-nonpositive-silently-gr | A15 | error by design | 0, -0.1, NaN, 2: CAMBValueError; 1e-30: CAMBError (segfaulted before 5b214ea) |
| pi-rk45-step-instability | C9 | fixed | 5_Horava_SolSyst_2 default/boost2 7.3e-4 = GR level |
| pk-extrap-slope-bias | B19 | fixed | identical to upstream |
| polterddot-nonflat-gr, t3-polterddot-curvature-gr | B9 | fixed | EFTflag = 0 non-flat vs CAMB <= 3.6e-15 |
| qs-mu-nan-vanishing-denominator | A23 (F-103) | fixed | no NaN, GR value 1 |
| qs-sigma-copy-of-qs-mu | B5 | fixed | QS_sigma = QS_sigma_gen |
| read-parameters-roundtrip-bytes | A8 | fixed | str/float returned, round trip unchanged |
| rgr-threshold-dimensionful-drops-de-perturbations | C3 | fixed | threshold 1e-7: sigma8 change 0 / 4.9e-9 |
| rph-kinetic-early-roundoff-ghost | C8 | fixed | healthy RPH model stable at stability time 1e-10 |
| rph-param-index-tensor-ode | A23 | fixed | last name/value correct |
| scg-adotoa-drops-massive-nu | B15 | fixed | adotoa vs CAMB H 4e-16 (second column of the script is obsolete after B11 changed grhoa2) |
| scg-higher-derivatives-wrong | B8 | fixed | all rows at the FD error (<= 1e-4) |
| set-dark-energy-ignored-with-active-eft | C11 | error by design | CAMBValueError "set_dark_energy cannot be used with an active EFTCAMB model" |
| set-params-cp-drops-eft-model, t3-set-params-cp-resets-eft-model | A6 | fixed | model kept (TT identical to control; EFTflag 3 after cp reuse) |
| stability-cache-cdotdot-uninitialised, uninit-cache-fields-garbage | A12 | fixed | verdicts and outputs independent of prior cache content |
| stability-gpinu-tot-tripled | A13 | fixed | ratio 1.0000 |
| stability-sampler-stale-y, stability-third-pass-stale-y | B4 | fixed | ghost window near a = 1 detected; pass 3 reaches aend - 1e-10 |
| t10a-forced-inttolboost-gr-limit | C10 | by design (documented) | EFT IntTolBoost >= 2 forced: LCDM-equivalent EFT vs EFTflag = 0 tensor TT 7.6e-2 at default, 1e-11 at matched IntTolBoost |
| t10a-upstream-fluid-w-near-minus1-rk45-failure | none | upstream | fails identically in B0 |
| t10a-upstream-get-background-windows-segfault | D (F-099) | upstream, still reproduces | WantCls = F with windows: segfault in HEAD and B0 |
| t3-3p1-success-unset | A4 | fixed (solver); refused by physics checks | the background solver works; audit ini and GR point refused ("theory unstable"), forced GR point within 4.3e-5 of EFTflag = 0 (open question) |
| t3-bg-table-resolution-gr-floor | B21 | fixed | 10x tables are now the default (identical rows); the remaining 3-4e-5 TT vs GR is the C9 floor |
| t3-counts-rsd-gr-branch-missing-de | B10 | fixed | counts C_l vs CAMB 1.4e-7 |
| t3-eft-curvature-ignored | C1 | error by design | omk != 0 runs refused |
| t3-eft-turnon-before-tauminn-lmax | A16 | fixed | l_max_scalar 2000 runs |
| t3-extgal-gr-point-background | B15 | fixed | S = 0: status ok, CAMB H = LCDM to 1e-16; the remaining 2.7e-5 in the adotoa row is the script's linear z interpolation of the reference H |
| t3-fr-rgr-switch-gr-approach | C3 | fixed | P(k) deviations smooth in B0 (RGR a 0.075-0.091) |
| t3-gamma1-mass-instability-near-gr | C4 | fixed | constant Gamma1 1e-2 ... 1e-9 accepted |
| t3-get-evolution-transfer-only-segfault | A3 | fixed | get_scale_evolution works after a transfer-only run |
| t3-pi-ic-singular-denominator | C2 | fixed | no spike; deviation -> 4e-6 as RPHmassP0 -> 0 |
| t3-stability-dtauda-without-matter | B11 | fixed (dtauda); still refused | no "dtauda is NaN"; OL-gamma with OLOmega0 <= 1e-6 and its GR point now refused as ghost (kinetic -5e-3 / -9.4e-3 at a = 1e-10; open question) |
| theta-h0-range-eft | D (F-139) | fixed | thetastar solved in both cases |
| theta-solve-aborts-unstable-bracket-end | A17 | fixed | H0 recovered to <= 4e-6 |
| wde-integral-unsupported-parametrizations | A9/B13 (F-047) | fixed | EFTwDE 6-11 run, <= 3.7e-11 from the constant reference |
| write-ini-drops-eftcamb | A7 | fixed | 47 EFTCAMB keys written, re-read TT identical |

Totals (97 scripts): 77 fixed (2 of them run by hand with a checked build of HEAD; 4 with a residual open question:
ic1-pi-attractor-pole, t3-3p1-success-unset, t3-stability-dtauda-without-matter, hdsk-h0-not-propagated theta),
13 error or behaviour by design, 2 partial (C9), 3 not applicable (deleted instrumented trees; items checked
otherwise), 2 upstream (also in plain CAMB). One crash found by this run was fixed (5b214ea).
