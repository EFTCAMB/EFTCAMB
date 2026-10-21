# EFTCAMB audit: final regression (HEAD vs pre-fix baseline)

HEAD = `eftcamb-audit` at 4819c8f (all approved fixes, A9 included). Baseline A0 = the pre-fix code (f6f4ba9) as frozen
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
- Common floor: most models with active EFT perturbations move by TT ~2e-5, EE ~1.1e-5, P(k) <= 1e-6. This is B3
  (pi''-dependent Einstein factors recomputed after the pi equation): the core-group tree with only B3 reverted
  (`fixA/attr/noB3`) reproduces A0 bit-for-bit on all 29 models tried and on the six 007p5/6/7 audit inis.
- Runtime: median HEAD/A0 0.98; total 2198 s -> 2268 s (+3%). Larger: 4_Des_F_R_LCDM_2 39.5 -> 109.5 s (C3, earlier
  turn-on 0.093 -> 0.016 where the f(R) pi equation is stiff), 008p5-beyond-horndeski(_2) 17.9/12.9 -> 28.3/19.1 s
  (B16, background from the model's own solution), 5_scg 20.2 -> 25.2 s (B15, massive nu in the ODE),
  007p5/007p6 +12-18%, 008p7-ext-galileon(_2) 204/177 s (failed after 58/31 s in A0).

## Changed models and responsible items

| model(s) | biggest change HEAD vs A0 | items |
|---|---|---|
| 2_PEFT_Omega_const_1/2, gamma3_*, gamma5_const_1, gamma5_exp_1, gamma4_const_1 | TT 1.2e-2 (Omega_const_1, gamma5_const_1), 5-8e-3 (gamma3), P(k) <= 7e-5, sigma8 <= 2.5e-5 | B3 |
| other 2_PEFT_* with EFT on, 4_Des_mc5e_*, 4_Des_F_R_LCDM_1, wCDM_1/2, 007p7-shiftsym(_2) | TT ~2e-5, EE 1.1e-5 (B3 floor); f(R) P(k) 2-4e-6 | B3; f(R) also stab group (C3 RGR, <= 3e-6) and B6 (Hdotdotdot, get_eft_functions only) |
| 3_RPH_MB_const_1 | P(k) x100 at a single k (the IC pole), TT 1.0, sigma8 5.7e-6 | C2 (pi IC regularisation; the A0 spike and zeros across k_c are gone, default and boost2 agree); B3 TT 7e-2 |
| 3_RPH_MB_OmegaDE_1, MK_const_1, KT/MK_OmegaDE_1, MK_power_1, KT_power_1 | TT 2.7e-3 / 6.5e-3 / 1e-4 | B3; stab group (C2/C8) <= 4.5e-5 |
| 4_Des_F_R_LCDM_2 | P(k) 0.13 (k >= 1/Mpc), sigma8 2.3e-4, PP 3.7e-3 | C3 (dimensionless, observable-based return to GR: turn-on 0.0931 -> 0.0162, converged answer); runtime x2.8 |
| 5_ADE_1 | 8.6e-5 C_l, P(k) 7.3e-5, H 3.2e-5 | B21 (10x background tables), C3 (turn-on 1e-10 -> 1.16e-6), B3 |
| 5_Horava_* | TT 0.49 (SolSyst_3), 0.38 (eta_4), 4.3e-2 (SolSyst_2); sigma8 <= 5.7e-3 | B7 (Hubble derivatives from the model's H(a), arXiv:1508.01787) and B3 |
| 5_Kmimic_1/2/3 | P(k) 18 / 530 / 0.86 (single-k IC poles removed), TT 0.94 / 0.92 / 0.22, sigma8 0.776 -> 0.616 etc. (2.1e-1) | C2 (IC regularisation, k-independent IC time), B1 (EFTpiDfunction), B3, B21 (10000-point table, c6b0fb0) |
| 5_Kmouflage_1/2/3 | TT 1.3e-2 / 7.7e-2 / 4e-4; H 7-8e-5 | B3; B21 (tables: 2.4e-4 C_l, 7e-5 H) |
| 5_hdsk_TG_1/2/3 | sigma8 0.14 / 0.14 / 0.022; C_l <= 2e-4 | C7 (EFT background H0 adopted); TG_2 also C2 (P(k) 7e-4, now agrees with TG_1 delta_phi) and B3 |
| 5_quint_* | C_l 1-2.4e-4, H 2-6e-5; quint_cubic P(k) 4.8e-4 | B21 (tables); C3 (turn-on 1e-8 -> 1-2.5e-3, frozen field, <= 3e-5); quint_cubic P(k) from the stab group |
| 5_scg | TT 0.21, P(k) 1.4, sigma8 7.2e-2, H 1e-2 | B15 + B8 + B23 (massive nu in the background, third derivatives, _dl literal) |
| 007p5-OL-gamma(_2), 007p6-OL-alpha(_2) (audit inis) | TT 1.5e-2 / 4.6e-3 / 4.9e-3 / 3.4e-3, P(k) <= 5.4e-4 | B3 (noB3 tree = A0 bit-for-bit) |
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
| 4_Des_F_R_LCDM_2 | ok | 3.2e-05 | 8.7e-06 | 5.1e-06 | 1.1e-03 | 3.7e-03 | 1.3e-01 | 2.3e-04 | 0 | 39.5 / 109.5 | bit | 1.9e-04 / 1.3e-01 / 0 |
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
