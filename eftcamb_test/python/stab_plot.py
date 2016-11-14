import CAMB_plots_lib.plot_stability_space as splt

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.tri as tri
import math

path = "../results/Stability_Results/"

EFT_1D_root                = path+"1_PureEFT_1P/"
PureEFT_wCDM_2D_root       = path+"2_PureEFT_wCDM_2D/"
PureEFT_PowerLaw_root      = path+"3_PureEFT_PowerLaw/"
PureEFT_Exponential_root   = path+"4_PureEFT_Exponential/"
Designer_1D_root           = path+"5_Designer_1D/"
Designer_2D_root           = path+"6_Designer_2D/"
Designer_RPH_1D_root       = path+"7_RPH_1P/"
Designer_RPH_2D_root       = path+"8_RPH_2P/"
Designer_RPH_2D_const_root = path+"9_RPH_2P_const/"
PureEFT_1D_Horndeski_root  = path+"10_PEFT_Horn_1D/"
Horava_1D_root             = path+"11_Horava_1D/"
Horava_2D_root             = path+"12_Horava_2D/"

outdir = "../results/Stability_Plots/"

""" Do 1 parameter EFT models """
         
models1D =[
         [EFT_1D_root+"1_Omega_constant.dat","$\Omega_0$","constant EFT $\Omega$", "Omega_constant"],
         [EFT_1D_root+"1_Omega_linear.dat","$\Omega_0$","linear EFT $\Omega$","Omega_linear"],
         [EFT_1D_root+"1_Gamma1_constant.dat","$\gamma_1^0$","constant EFT $\gamma_1$","Gamma1_constant"],
         [EFT_1D_root+"1_Gamma1_linear.dat","$\gamma_1^0$","linear EFT $\gamma_1$","Gamma1_linear"],
         [EFT_1D_root+"1_Gamma2_constant.dat","$\gamma_2^0$","constant EFT $\gamma_2$","Gamma2_constant"],
         [EFT_1D_root+"1_Gamma2_linear.dat","$\gamma_2^0$","linear EFT $\gamma_2$","Gamma2_linear"],
         [EFT_1D_root+"1_Gamma3_constant.dat","$\gamma_3^0$","constant EFT $\gamma_3$","Gamma3_constant"],
         [EFT_1D_root+"1_Gamma3_linear.dat","$\gamma_3^0$","linear EFT $\gamma_3$","Gamma3_linear"],
         [EFT_1D_root+"1_Gamma4_constant.dat","$\gamma_4^0$","constant EFT $\gamma_4$","Gamma4_constant"],
         [EFT_1D_root+"1_Gamma4_linear.dat","$\gamma_4^0$","linear EFT $\gamma_4$","Gamma4_linear"],
         [EFT_1D_root+"1_Gamma5_constant.dat","$\gamma_5^0$","constant EFT $\gamma_5$","Gamma5_constant"],
         [EFT_1D_root+"1_Gamma5_linear.dat","$\gamma_5^0$","linear EFT $\gamma_5$","Gamma5_linear"],
         [EFT_1D_root+"1_Gamma6_constant.dat","$\gamma_6^0$","constant EFT $\gamma_6$","Gamma6_constant"],
         [EFT_1D_root+"1_Gamma6_linear.dat","$\gamma_6^0$","linear EFT $\gamma_6$","Gamma6_linear"],  
         [Designer_1D_root+"1_des_5e_wCDM.dat","$w_0$","minimally coupled quintessence on wCDM","des_5e_wCDM"],
         [Designer_1D_root+"1_des_fR_LCDM.dat","$B_0$","designer f(R) on $\Lambda$CDM","des_fR_LCDM"],
         [Designer_RPH_1D_root+"1_RPH_Braiding_const.dat","$a_B^0$","RPH braiding constant","RPH_Braid_const"],
         [Designer_RPH_1D_root+"1_RPH_Kineticity_const.dat","$a_K^0$","RPH kineticity constant","RPH_Kinet_const"],
         [Designer_RPH_1D_root+"1_RPH_MassP_const.dat","$M^0$","RPH Planck Mass constant","RPH_MassP_const"],
         [Designer_RPH_1D_root+"1_RPH_Tensor_const.dat","$a_T^0$","RPH tensor constant","RPH_Tensor_const"],   
         [PureEFT_1D_Horndeski_root+"1_Gamma3_Hor_constant.dat", "$\gamma_3^0$","PureEFT Horndeski $\gamma_3$ const","Gamma3_Hor_const"],
         [PureEFT_1D_Horndeski_root+"1_Gamma3_Hor_linear.dat", "$\gamma_3^0$","PureEFT Horndeski $\gamma_3$ linear","Gamma3_Hor_linear"],
         [Horava_1D_root+"1_Horava_xi.dat", "$xi$","Horava xi","Horava_xi"],
         [Horava_1D_root+"2_Horava_lambda.dat", "$\lambda$","Horava lambda","Horava_lambda"],
         [Horava_1D_root+"3_Horava_eta.dat", "$\eta$","Horava eta","Horava_eta"],
         ]

for M in models1D:
    temp = splt.stability_1D(M[0],M[1],M[2])
    temp.stab_plot()
    plt.savefig(outdir+M[3],aspect='auto',bbox_inches='tight')
    plt.clf()

""" Do 2 parameters EFT models """

models2D = [
          [PureEFT_wCDM_2D_root+"1_Omega_const_wCDM.dat", "$\Omega_0$", "$w_0$", "constant EFT $\Omega$ on wCDM", "Omega_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Omega_linear_wCDM.dat", "$\Omega_0$", "$w_0$", "linear EFT $\Omega$ on wCDM", "Omega_linear_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma1_const_wCDM.dat", "$\gamma_1^0$", "$w_0$", "constant EFT $\gamma_1$ on wCDM", "Gamma1_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma1_linear_wCDM.dat", "$\gamma_1^0$", "$w_0$", "linear EFT $\gamma_1$ on wCDM", "Gamma1_linear_wCDM"],  
          [PureEFT_wCDM_2D_root+"1_Gamma2_const_wCDM.dat", "$\gamma_2^0$", "$w_0$", "constant EFT $\gamma_2$ on wCDM", "Gamma2_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma2_linear_wCDM.dat", "$\gamma_2^0$", "$w_0$", "linear EFT $\gamma_2$ on wCDM", "Gamma2_linear_wCDM"], 
          [PureEFT_wCDM_2D_root+"1_Gamma3_const_wCDM.dat", "$\gamma_3^0$", "$w_0$", "constant EFT $\gamma_3$ on wCDM", "Gamma3_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma3_linear_wCDM.dat", "$\gamma_3^0$", "$w_0$", "linear EFT $\gamma_3$ on wCDM", "Gamma3_linear_wCDM"], 
          [PureEFT_wCDM_2D_root+"1_Gamma4_const_wCDM.dat", "$\gamma_4^0$", "$w_0$", "constant EFT $\gamma_4$ on wCDM", "Gamma4_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma4_linear_wCDM.dat", "$\gamma_4^0$", "$w_0$", "linear EFT $\gamma_4$ on wCDM", "Gamma4_linear_wCDM"], 
          [PureEFT_wCDM_2D_root+"1_Gamma5_const_wCDM.dat", "$\gamma_5^0$", "$w_0$", "constant EFT $\gamma_5$ on wCDM", "Gamma5_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma5_linear_wCDM.dat", "$\gamma_5^0$", "$w_0$", "linear EFT $\gamma_5$ on wCDM", "Gamma5_linear_wCDM"], 
          [PureEFT_wCDM_2D_root+"1_Gamma6_const_wCDM.dat", "$\gamma_6^0$", "$w_0$", "constant EFT $\gamma_6$ on wCDM", "Gamma6_constant_wCDM"],
          [PureEFT_wCDM_2D_root+"1_Gamma6_linear_wCDM.dat", "$\gamma_6^0$", "$w_0$", "linear EFT $\gamma_6$ on wCDM", "Gamma6_linear_wCDM"],
          
          [PureEFT_PowerLaw_root+"1_Omega_PL.dat", "$\Omega_0$", "$s$", "power law EFT $\Omega$", "Omega_pow"],
          [PureEFT_PowerLaw_root+"1_Gamma1_PL.dat", "$\gamma_1^0$", "$s$", "power law EFT $\gamma_1$", "Gamma1_pow"],
          [PureEFT_PowerLaw_root+"1_Gamma2_PL.dat", "$\gamma_2^0$", "$s$", "power law EFT $\gamma_2$", "Gamma2_pow"],
          [PureEFT_PowerLaw_root+"1_Gamma3_PL.dat", "$\gamma_3^0$", "$s$", "power law EFT $\gamma_3$", "Gamma3_pow"],
          [PureEFT_PowerLaw_root+"1_Gamma4_PL.dat", "$\gamma_4^0$", "$s$", "power law EFT $\gamma_4$", "Gamma4_pow"],
          [PureEFT_PowerLaw_root+"1_Gamma5_PL.dat", "$\gamma_5^0$", "$s$", "power law EFT $\gamma_5$", "Gamma5_pow"],
          [PureEFT_PowerLaw_root+"1_Gamma6_PL.dat", "$\gamma_6^0$", "$s$", "power law EFT $\gamma_6$", "Gamma6_pow"],
          
          [PureEFT_Exponential_root+"1_Omega_exp.dat", "$\Omega_0$", "$s$", "exponential EFT $\Omega$", "Omega_exp"],
          [PureEFT_Exponential_root+"1_Gamma1_exp.dat", "$\gamma_1^0$", "$s$", "exponential EFT $\gamma_1$", "Gamma1_exp"],
          [PureEFT_Exponential_root+"1_Gamma2_exp.dat", "$\gamma_2^0$", "$s$", "exponential EFT $\gamma_2$", "Gamma2_exp"],
          [PureEFT_Exponential_root+"1_Gamma3_exp.dat", "$\gamma_3^0$", "$s$", "exponential EFT $\gamma_3$", "Gamma3_exp"],
          [PureEFT_Exponential_root+"1_Gamma4_exp.dat", "$\gamma_4^0$", "$s$", "exponential EFT $\gamma_4$", "Gamma4_exp"],
          [PureEFT_Exponential_root+"1_Gamma5_exp.dat", "$\gamma_5^0$", "$s$", "exponential EFT $\gamma_5$", "Gamma5_exp"],
          [PureEFT_Exponential_root+"1_Gamma6_exp.dat", "$\gamma_6^0$", "$s$", "exponential EFT $\gamma_6$", "Gamma6_exp"],
          
          [Designer_2D_root+"1_des_5e_CPL.dat", "$w_0$", "$w_a$", "minimally coupled quintessence on CPL","des_5e_CPL"],
          [Designer_2D_root+"1_des_fR_wCDM.dat", "$B_0$", "$w_0$", "designer f(R) on wCDM","des_fR_wCDM"],
          
          [Designer_RPH_2D_root+"1_Braiding_PL.dat", "$a_B^0$", "$s$", "RPH braiding power law","RPH_Braid_PL"],
          [Designer_RPH_2D_root+"1_Kineticity_PL.dat", "$a_K^0$", "$s$", "RPH kineticity power law","RPH_Kinet_PL"],
          [Designer_RPH_2D_root+"1_MassP_PL.dat", "$M^0$", "$s$", "RPH Planck Mass power law","RPH_MassP_PL"],
          [Designer_RPH_2D_root+"1_Tensor_PL.dat", "$a_T^0$", "$s$", "RPH tensor power law","RPH_Tensor_PL"],
          
          [Designer_RPH_2D_const_root+"1_Kinet_Braiding.dat", "$a_K^0$", "$a_B^0$", "RPH constant kinet and braiding","RPH_Kinet_Braid"],
          [Designer_RPH_2D_const_root+"1_Kinet_Tensor.dat", "$a_K^0$", "$a_T^0$", "RPH constant kinet and tensor","RPH_Kinet_Tensor"],
          [Designer_RPH_2D_const_root+"1_MassP_Braiding.dat", "$M_0$", "$a_B^0$", "RPH Mass and Braiding","RPH_Mass_Braid"],
          [Designer_RPH_2D_const_root+"1_MassP_Kinet.dat", "$M_0$", "$a_K^0$", "RPH Mass and Kinet","RPH_MassP_Kinet"],
          [Designer_RPH_2D_const_root+"1_MassP_Tensor.dat", "$M_0$", "$a_T^0$", "RPH Mass and Tensor","RPH_Mass_Tensor"],
          [Designer_RPH_2D_const_root+"1_Tensor_Braiding.dat", "$a_T^0$", "$a_B^0$", "RPH Tensor and Braiding","RPH_Tensor_Braiding"],
          
          [Horava_2D_root+"1_Horava_xi_lambda.dat", "$xi$", "$\lambda$", "Horava xi lambda","Horava_xi_lambda"],
          [Horava_2D_root+"2_Horava_lambda_eta.dat", "$\lambda$", "$\eta$", "Horava lambda eta","Horava_lambda_eta"],
          [Horava_2D_root+"3_Horava_xi_eta.dat", "$xi$", "$\eta$", "Horava xi eta","Horava_xi_eta"],
          [Horava_2D_root+"4_Horava_xi_lambda_SolSyst.dat", "$xi$", "$\lambda$", "Horava xi lambda with Sol Syst","Horava_xi_lambda_SolSyst"],
            
          ]


for M in models2D:
    temp = splt.stability_2D(M[0],[M[1],M[2]],M[3])
    temp.stab_plot()
    plt.savefig(outdir+M[4],aspect='auto',bbox_inches='tight')
    plt.clf()








