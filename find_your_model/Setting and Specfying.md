# Setting Parameters in the EFT Structure

#### There are two methods to set parameters in the EFT structure:

### 1. **Case-by-Case Definition**

   <span style="font-size:20px;"> This method is suitable for specific models and functions in the structure. Parameters have physical significance within their corresponding models.

   <span style="font-size:20px;"> Exact examples will also be provided in the examples folder, so we only list their names here. For several specific models, we also list their corrresponding reference papers. *Please add the citation when you use these models.*

   * <span style="font-size:19px;"> Flag for different $w_{DE}$ parametrizations (LCDM, wCDM, CPL, JBL, Taylor or TurningPoint) in PureEFTmodel (**EFTflag**=1): **EFTwDE**. </span>

        - **EFTwDE** = 0&emsp; ->&emsp; $w_{DE} = -1$
        - **EFTwDE** = 1&emsp; ->&emsp; $w_{DE} = w_0$
        - **EFTwDE** = 2&emsp; ->&emsp; $w_{DE} = w_0 + w_a(1-a)$
        - **EFTwDE** = 3&emsp; ->&emsp; $w_{DE} = w_0 + w_a (1-a) a^{(n-1)}$
        - **EFTwDE** = 4&emsp; ->&emsp; $w_{DE} = w_0 + w_a (a_t-a)^2$
        - **EFTwDE** = 5&emsp; ->&emsp; $w_{DE} = w_0 + w_a a + \frac12 w_2 a^2 + \frac16 w_3 a^3$
        - **EFTwDE** equals from 6 to 11 share the same logic of specifying functions decribed below, see **Specifying Functions** for more details.


      <span style="font-size:19px;">The parameters above can be fixed with the flags:</span>
       - $w_0$->**EFTw0** ,  $w_a$->**EFTwa** , $n$->**EFTwn**, $a_t$->**EFTwat** , $w_2$->**EFTw2** , $w_3$->**EFTw3**

      The flags for parametrizations of PureEFTmodel functions:
      - **PureEFTmodelOmega**, **PureEFTmodelGamma1**, **PureEFTmodelGamma2**,**PureEFTmodelGamma3**, **PureEFTmodelGamma4**, **PureEFTmodelGamma5**, **PureEFTmodelGamma6**.

      And their function names:
      - **EFTOmega**,**EFTGamma1**,**EFTGamma2**,**EFTGamma3**,**EFTGamma4**,**EFTGamma5**,**EFTGamma6**

      ---


   * <span style="font-size:19px;">Flags and Parameter names of **AltParEFTmodel** (**EFTflag**=2):</span>
      <br>
      * **AltParEFTmodel** = 1:
         - Flags for different $w_{DE}$ parametrizations: **RPHwDE**.
         So the parameter names are the same with **EFTwDE**, except changing ***"EFT"*** instead of ***"RPH"***.
         - Flags for parametrizations of RPH functions:
         **RPHalphaMmodel/RPHmassPmodel, RPHkineticitymodel, RPHbraidingmodel, RPHtensormodel**.
         ***(Parametrizations name rule see the text below.)***
         Their function names are **RPHalphaM, RPHmassP, RPHkineticity, RPHbraiding, RPHtenso** and Latex format are $\alpha^{\rm M},\tilde{M},\alpha^{\rm K},{\alpha^{\rm B}},{\alpha^{\rm T}}$, respectively.
         - With **RPHusealphaM** = T the Planck mass is integrated from $\alpha^{\rm M}$ and **RPH_M0** sets
           $1 + M_0 = M_{\rm eff}^2/m_0^2$ where the integration starts: today if **RPHintegratefromtoday** = T, at the
           start of the background (**model_background_a_ini**) otherwise. In the second case today's value differs
           (constant $\alpha^{\rm M} = 0.1$ with $M_0 = 0.1$ gives $M_{\rm eff}^2(a=1) = 6.94$).
      <br>
      * **AltParEFTmodel** = 2
         - Flags: **OLLambdamodel, OLOmegamodel, OLGamma1model, OLGamma2model, OLGamma3model**
         - Names: **OLLambda, OLOmega, OLGamma1, OLGamma2, OLGamma3**
         - Latex format: $\tilde{\lambda}, \Omega, \gamma_1, \gamma_2, \gamma_3$
         - Definitions: $\Lambda = \Lambda_0 (1 + \tilde{\lambda}(a))$ with $\Lambda_0 = -3 H_0^2 (1 - \Omega_m)$, the
           cosmological constant of $\Lambda$CDM with the same matter; $\gamma_4 = -\gamma_3$, $\gamma_5 = \gamma_3/2$
           (Horndeski); the background follows from $H(a=1) = H_0$. Note: $\Omega \neq 0$ with $\tilde{\lambda} = 0$ is
           not $\Lambda$CDM, flatness generates a dust-like $c = \frac{3}{2} H_0^2 \Omega_0 a^{-3}$ (a ghost for
           $\Omega_0 < 0$).
         - Reference papers: arXiv:1703.05297 (App. A), arXiv:2407.15832.
      <br>
      * **AltParEFTmodel** = 3
        - Flags: **OLLambdamodel, OLmassPmodel, OLkineticitymodel, OLbraidingmodel, OLtensormodel**
        - Names: **OLLambda, OLmass, OLkineticity, OLbraiding, OLtensor**
        - Latex format: $\tilde{\lambda}, \tilde{M}, \alpha^{\rm K}, \alpha^{\rm B}, \alpha^{\rm T}$
        - Definitions: $\tilde{\lambda}$ as for **AltParEFTmodel** = 2; $\tilde{M} = M_*^2/m_0^2 - 1$; $\alpha^{\rm B}$ is in
          the EFTCAMB convention, $-1/2$ times the $\alpha_B$ of Bellini-Sawicki and hi_class (enter $-\alpha_B/2$). As
          for **AltParEFTmodel** = 2, a constant $\tilde{M} \neq 0$ with $\tilde{\lambda} = 0$ is not $\Lambda$CDM:
          flatness generates a dust-like $c = \frac{3}{2} H_0^2 \tilde{M}_0 a^{-3}$.
        - Reference papers: arXiv:1703.05297 (App. A, background), arXiv:1405.3590 (alpha convention).
      <br>


      * **AltParEFTmodel** = 4
         Flag for different $w_{DE}$ parametrizations: **EFTwDE**.
         - Names: **Shift_Symmetric_alphaB0, Shift_Symmetric_m, Shift_Symmetric_alphaK0**
         - Latex format: $\hat{\alpha}_{B}, m, \alpha_{K0}$
         - Conventions: $\alpha_B = \hat{\alpha}_B (H_0/H)^{4/m}$ with $\hat{\alpha}_B$ the Bellini-Sawicki amplitude of
           arXiv:2103.11195, i.e. $-2$ times the $\alpha_{B,0}$ of the reference paper below (its best fit
           $\alpha_{B,0} = -0.30$ is **Shift_Symmetric_alphaB0** = 0.60); $m > 0$ is required unless $\hat{\alpha}_B = 0$;
           $\alpha_K = \alpha_{K0}\, a$, not constant as in both papers.
         **Reference Paper:**
         - [1] Albuquerque, I. S., Frusciante, N., Pace, F., and Schimd, C., “Spherical collapse and halo abundance in shift-symmetric Galileon theory”, <i>Physical Review D</i>, vol. 109, no. 2, Art. no. 023535, APS, 2024. [doi:10.1103/PhysRevD.109.023535](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.023535).

     ---

   * <span style="font-size:19px;"> Flags and Parameter names of designer mapping EFT (**EFTflag**=3): **MappingEFTmodel** </span>
      * **MappingEFTmodel** = 1: f(R)
         Needed parameters:
         - Names: **EFTB0**
         - Latex format: $\rm B_0$

         The parametrizations of wDE is similar with PureEFT case, with the prefix ***"EFT"***.

      * **MappingEFTmodel** = 3: Freezing Gravity

         <!-- The parametrizations of wDE is similar with PureEFT case, with the prefix ***"EFT"***. -->

         **Reference Paper:**
         - [1]Yao, Z., Ye, G., and Silvestri, A., “A General Model for Dark Energy Crossing the Phantom Divide”, <i>arXiv e-prints</i>, Art. no. arXiv:2508.01378, 2025. [doi:10.48550/arXiv.2508.01378](https://arxiv.org/abs/2508.01378).

     ---

   * <span style="font-size:19px;">Flags and Parameter names of full EFT mapping (**EFTflag**=4): **FullMappingEFTmodel** </span>
       <br>
      * **FullMappingEFTmodel** = 1: Horava
         Needed parameters:
         - Names: **Horava_eta, Horava_xi, Horava_lambda** (if **HoravaSolarSystem** = 1, then **Horava_xi** is not needed)
         - Latex format: $\eta_{\rm Ho\v rava}, \xi_{\rm Ho\v rava} - 1, \lambda_{\rm Ho\v rava} - 1$
         - Convention: deviations from GR (all zero is GR) in the notation of [1]: **Horava_eta** = $\eta$,
           **Horava_xi** = $\xi - 1$, **Horava_lambda** = $\lambda - 1$. The example $\lambda = 1.4$, $\xi = 0.9$ of [1] is
           **Horava_lambda** = 0.4, **Horava_xi** = -0.1 (entered as is it runs $\lambda = 2.4$, $\xi = 1.9$). The solar-system
           class is **Horava_xi** = **Horava_eta**/2, i.e. $\eta = 2\xi - 2$.

         Reference Paper:
         - [1] Frusciante, N., Raveri, M., Vernieri, D., Hu, B., and Silvestri, A., “Hořava Gravity in the Effective Field Theory formalism: From cosmology to observational constraints”, <i>Physics of the Dark Universe</i>, vol. 13, pp. 7–24, 2016. [doi:10.1016/j.dark.2016.03.002](https://www.sciencedirect.com/science/article/pii/S2212686416300139?via%3Dihub).
         - [2] Frusciante, N. and Benetti, M., “Cosmological constraints on Hořava gravity revised in light of GW170817 and GRB170817A and the degeneracy with massive neutrinos”, <i>Physical Review D</i>, vol. 103, no. 10, Art. no. 104060, APS, 2021. [doi:10.1103/PhysRevD.103.104060](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.104060).
       <br>
      * **FullMappingEFTmodel** = 2: Acoustic Dark Energy (ADE)
         Needed parameters:
         - Names: **cs2, Log_ac, f_ac, p, wf**
         - Latex format: $c_s^2, \log(a_c), f(a_c), p, w_f$
         - Model of arXiv:1905.12618: $1 + w_{\rm ADE} = (1 + w_f)/[1 + (a_c/a)^{3(1+w_f)/p}]^p$, $f(a_c)$ the ADE fraction
           of the total density at $a_c$, $c_s^2$ its rest-frame sound speed. The default **p** = 1 differs from the
           paper's $p = 1/2$ (the test model `5_ADE_1` sets 0.5).
       <br>
      * **FullMappingEFTmodel** = 3: K-mouflage
         Needed parameters:
         - Names: **alphaU, gammaU, m, eps2_0, gammaA**
         - Latex format: $\alpha_{U}, \gamma_{U}, m, \epsilon_{2,0}, \gamma_{A}$

         Reference Paper:
         - [1] Benevento, G., “K-mouflage imprints on cosmological observables and data constraints”, <i>Journal of Cosmology and Astroparticle Physics</i>, vol. 2019, no. 5, Art. no. 027, IOP, 2019. [doi:10.1088/1475-7516/2019/05/027.](https://iopscience.iop.org/article/10.1088/1475-7516/2019/05/027)
       <br>
      * **FullMappingEFTmodel** = 4: Quintessence
         There are two flags within Quintessence model: **potential_model** and **drag_initial_conditions**.
         - **potential_model** select the parametrization model for the quintessence potential.
         - Conventions: the field is in units of the reduced Planck mass and the potential in units of
           $3H_0^2(1-\Omega_c-\Omega_b-\Omega_\nu)/(8\pi G)$, which does not subtract photons and massless neutrinos.
           The $\Lambda$CDM point (a constant potential, **potential_model** = 1) is therefore
           $V_0 = \Omega_\Lambda/(1-\Omega_c-\Omega_b-\Omega_\nu)$ (e.g. 0.999895), not 1.
         - **phidot_ini** is $d\phi/d\ln a$ at the initial time (with **drag_initial_conditions** = F); with
           **drag_initial_conditions** = T the initial velocity has $d^2\phi/d\ln a^2 = 0$ (field dragged by the
           Hubble friction). The initial field is found by shooting on $H_0$ within [**field_min**, **field_max**].
       <br>
      * **FullMappingEFTmodel** = 5: Beyond Horndeski
         Needed parameters:
         - Names: **Beyond_Horndeski_x10, Beyond_Horndeski_x30, Beyond_Horndeski_x40**
         - Latex format: $x_1^0，x_3^0，x_4^0$

         Reference Paper:
         - [1] Peirone, S., Benevento, G., Frusciante, N., and Tsujikawa, S., “Cosmological constraints and phenomenology of a beyond-Horndeski model”, <i>Physical Review D</i>, vol. 100, no. 6, Art. no. 063509, APS, 2019. [doi:10.1103/PhysRevD.100.063509](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.063509).
       <br>
      * **FullMappingEFTmodel** = 6: Scaling Cubic Galileon

         Needed parameters:
         - Names: **Scaling_Cubic_A, Scaling_Cubic_beta1, Scaling_Cubic_beta2, Scaling_Cubic_lambda**
         - Latex format: $A, \beta_1, \beta_2, \lambda$

         Note: the positivity bounds (**EFT_positivity_bounds**) and the Minkowski-limit check (**EFT_minkowski_limit**)
         do not apply to this model. They test the Horndeski theory reconstructed from the EFT functions with $G_i$
         analytic at $X = 0$, while here $G_3 = A \ln(X e^{\lambda\phi})$: with either check on, the published models
         and the test model `5_scg` are refused at every time. Keep both off.

         Reference Paper:
         - [1] Albuquerque, I. S., Frusciante, N., and Martinelli, M., “Constraining cosmological scaling solutions of a Galileon field”, <i>Physical Review D</i>, vol. 105, no. 4, Art. no. 044056, APS, 2022. [doi:10.1103/PhysRevD.105.044056](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.105.044056).
       <br>
      * **FullMappingEFTmodel** = 7: Extended Galileon
         Reference Paper:
         - [1]Frusciante, N., Peirone, S., Atayde, L., and De Felice, A., “Phenomenology of the generalized cubic covariant Galileon model and cosmological bounds”, <i>Physical Review D</i>, vol. 101, no. 6, Art. no. 064001, APS, 2020. [doi:10.1103/PhysRevD.101.064001.](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.064001).
       ---
   * <span style="font-size:20px;">Flags and Parameter names of Horndeski Module (**EFTflag**=5): </span>
      <span style="font-size:18px;">This module allows the user to work directly with any covariant Lagrangian belonging to the Horndeski class. You can map your own theory within Horndeski Gravity with the coeffcients. Here we show the usages of Early Dark Energy and Thawling Gravity as two examples in the *example_notebooks* file.</span>

      Reference Paper:
       - [1]Ye, G., Martinelli, M., Hu, B., and Silvestri, A., “Hints of Nonminimally Coupled Gravity in DESI 2024 Baryon Acoustic Oscillation Measurements”, <i>Physical Review Letters</i>, vol. 134, no. 18, Art. no. 181002, APS, 2025. [doi:10.1103/PhysRevLett.134.181002.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.181002).
       - [2]Ye, G., “Bridge the Cosmological Tensions with Thawing Gravity”, <i>arXiv e-prints</i>, Art. no. arXiv:2411.11743, 2024. [doi:10.48550/arXiv.2411.11743.](https://arxiv.org/abs/2411.11743)







### 2. **Default Definition**

   <span style="font-size:19px;"> If there is no specific definition for the model, a default name will be used. This name is constructed as "Function name" + index, where "Function name" is defined in the model flag, and the index ranges from 0 to \( N-1 \), where \( N \) is the number of parameters. The LaTeX format will be: **$\text{FunctionName}_{\rm index}$** . </span>

## Specifying functions

<span style="font-size:19px;"> In the flowchart, you will notice several functions labeled with *Specifying Functions*. We have implemented a number of parametrization methods within EFTCAMB. The functions can now be defined as both functions of the scale factor and the DE energy fraction $\Omega_{\rm DE}$. For Horndeski modual,  the Horndeski functions can also be defines as functions of scaler field $\phi$. The flag numbers and corresponding parameter names are listed here. Still, if there is no specific definition for a parameter name, the default definition will be used. </span>

<span style="font-size:19px;"> **We use $\Omega$ in pureEFT as an example** </span>

   Flags for specifying $\Omega$: **PureEFTmodelOmega**

   - Zero: **PureEFTmodelOmega** = 0&emsp; ->&emsp; $\Omega(a) = 0$
    <br>
   - Constant: **PureEFTmodelOmega** = 1&emsp; ->&emsp; $\Omega(a) = \Omega_0$
    <br>
   - Linear: **PureEFTmodelOmega** = 2&emsp; ->&emsp; $\Omega(a) = \Omega_0 a$
    <br>
   - Power Law: **PureEFTmodelOmega** = 3&emsp; ->&emsp; $\Omega(a) = \Omega_0 a^s$.
   Parameter Names: $\Omega_0$-> **EFTOmega0** , $s$ -> **EFTOmegaExp**.
   Latex format label: ($\Omega_0$, $\Omega_n$).
    <br>
   - Exponential: **PureEFTmodelOmega** = 4&emsp; ->&emsp; $\Omega(a) = \exp(\Omega_0 a^s) -1$.
   Parameter Names: $\Omega_0$-> **EFTOmega0** , $s$ -> **EFTOmegaExp**.
   Latex format label: ($\Omega_0$, $\Omega_n$).
    <br>
   - Taylor series: **PureEFTmodelOmega** = 5&emsp; ->&emsp; $\Omega(a) = \sum_{i=0}^{n} \Omega_{n}(a-a_0)^n$.
   **EFTOmega_Taylor_order** define the order and **EFTOmegaa0** define the expansion time.
   Names of the other parameters: Default defination.
    <br>
   - Pade series: **PureEFTmodelOmega** = 6&emsp; ->&emsp; $\Omega(a) = \frac{\sum_{i=0}^{n} \Omega_{\rm up,n}(a-a_0)^n}{1+\sum_{i=1}^{m} \Omega_{\rm down,m}(a-a_0)^m} $
   **EFTOmega_Pade_order_N** and **EFTOmega_Pade_order_D** define the order on the numerator and denominator. **EFTOmegaa0** define the expansion time.
   Names of other parameters: Default defination with FunctionName+N for numerator and FunctionName+D for denominator.
    <br>
   - Fourier: **PureEFTmodelOmega** = 7&emsp; ->&emsp; $\Omega(a) = a_0 + \sum_{i=1}^{n}cos[2\pi n(a-a_p)] + \sum_{i=1}^{m}sin[2\pi m(a-a_p)]$
   **EFTOmega_Fourier_order_cos** and **EFTOmega_Fourier_order_sin** define the order for cos and sin.
   **EFTOmega_a0** define the expansion time. **EFTOmega_phase** define the phase.
   Names of other parameters: Default defination with FunctionName+_a for numerator and FunctionName+_b for denominator.
   <br>
   - Steplog: **PureEFTmodelOmega** = 8&emsp; ->&emsp; $ \Omega(a) = \frac{(v_2-v_1)log\frac{a}{a_{\rm T}}}{2\delta\sqrt{1+(\frac{log\frac{a}{a_{\rm T}}}{\delta}})^2}+\frac{v_1+v_2}{2}$
   Names of the parameters: **EFTOmega_v1, EFTOmega_v2, EFTOmega_at, EFTOmega_delta**
   <br>
   - Spline (spline_rec): **PureEFTmodelOmega** = 9&emsp;
   Spline interpolation is a method used to construct a smooth curve that passes through a given set of data points. Here we implement the cubic and fifth-order splines for non-parametrizations of the functions you want to specify.
   This is the natural cubic spline of the nodes with approximate derivatives: each derivative is a separate cubic spline through derivatives estimated at the nodes, set to zero at the first and last node. It is not the derivative of the value spline. Use Spline3 (flag 12) for exact derivatives.
   **EFTOmega_Spline_Pixels** defines the number N of chosen pixels.
   Name of the parameters: **EFTOmegaxn,n$\in${1,N}** are the scale coordinates, **EFTOmegavn,n$\in${1,N}** are their values.
   **EFTOmega_null_value** is the value of the function outside the nodes (default 0); there the derivatives are 0.
   *Notes: [1] Spline x coordinate must be increasing. [2] N to be at least 4.*
   <br>
   - Spline5: **PureEFTmodelOmega** = 10&emsp;
   Fifth order spline.
   <br>
   - Exponential_Parametrization_2_1D: **PureEFTmodelOmega** = 11&emsp; ->&emsp; $\Omega(a) = \Omega_0e^{-sa}$.
   The other parametrization method of exponential.
   Names of the parameters: Default defination
   <br>
   - Spline3: **PureEFTmodelOmega** = 12&emsp;
   Natural cubic spline of the nodes with exact derivatives: the first, second and third derivatives are those of the value spline, the fourth derivative is zero.
   Same parameters as Spline (flag 9): **EFTOmega_Spline_Pixels**, **EFTOmegaxn**, **EFTOmegavn**, **EFTOmega_null_value**.
   *Notes: [1] Spline x coordinate must be increasing. [2] N to be at least 4.*
