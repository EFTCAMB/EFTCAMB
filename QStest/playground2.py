# initial imports:
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.interpolate import interp1d
from getdist import loadMCSamples
from getdist import MCSamples
from utilities.settings import *
camb_path = os.path.realpath(os.path.join(os.getcwd(), here+'/../'))
sys.path.insert(0, camb_path)
from camb import model
import camb

camb.set_feedback_level(20)

###############################################################################
# define parameters:

# hardcoded lensing bf
bestfit = {}
bestfit["ombh2"] = 0.2266705E-01
bestfit["omch2"] = 0.1169748E+00
bestfit["As"] = 0.2092396E+01/10**9
bestfit["H0"] = 0.6882996E+02
bestfit["tau"] = 0.5617335E-01
bestfit["ns"] = 0.9742312E+00
bestfit["mnu"] = 0.6000000E-01
bestfit["mnu"] = 0.06
bestfit["nnu"] = 3.046

#minimal const omega
# eftcamb_params = {'EFTflag': 1,
#                   'PureEFTmodel': 1,
#                   'PureEFTmodelOmega': 4,
#                   'EFTwDE': 0,
#                   'EFTOmega0': 1.05,
#                   'EFTOmegaExp': 2.05,
#                   'PureEFTmodelGamma1': 0,
#                   'PureEFTmodelGamma2': 0,
#                   'PureEFTmodelGamma3': 0,
#                   'PureEFTmodelGamma4': 0,
#                   'PureEFTmodelGamma5': 0,
#                   'PureEFTmodelGamma6': 0,
#                   'PureEFTHorndeski': True,
#                   "EFT_additional_priors": False,
#                   'feedback_level': 20,
#                   'EFTCAMB_turn_on_time': 1.e-12,
#                   'EFTCAMB_do_QSA': True
#                   }

#5_Kmimic_1
# eftcamb_params = {'EFTflag': 4,
#                   'EFTflag': 4,
#                   'FullMappingEFTmodel': 3,
#                   'Kmimic': True,
#                   'alphaU': 0.2,
#                   'gammaU': 1.,
#                   'm': 3.0,
#                   'eps2_0': 0.01,
#                   'gammaA': 5.0,
#                   'omnuh2': 0.0,
#                   'feedback_level': 20,
#                   'EFTCAMB_turn_on_time': 1.e-12,
#                   'EFTCAMB_do_QSA': True
#                   }
#

# 2_PEFT_gamma4_power_2
eftcamb_params = {  'EFTflag' : 1,
                    'PureEFTmodel'     : 1,
                    'PureEFTHorndeski' : True,
                    'EFTwDE'             : 0,
                    'PureEFTmodelOmega'  : 1,
                    'PureEFTmodelGamma1' : 1,
                    'PureEFTmodelGamma2' : 1,
                    'PureEFTmodelGamma3' : 0,
                    'PureEFTmodelGamma4' : 0,
                    'PureEFTmodelGamma5' : 0,
                    'PureEFTmodelGamma6' : 0,

                    'EFTGamma10' : 0.5,
                    'EFTGamma20' : 0.01,
                    #'EFTGamma30' : 0.001,
                    # 'EFTGamma40' : 0.001,
                    'EFTOmega0' : 0.05,
                    "EFT_ghost_stability" : False,
                    "EFT_gradient_stability" : False,
                  'feedback_level': 20,
                  'EFTCAMB_turn_on_time': 1.e-12,
                   'EFTCAMB_do_QSA': True
                  }


out_dir = results_dir+'/playground/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

#
pars_QSA = camb.set_params(**eftcamb_params, **bestfit, DoLateRadTruncation=True)
# test 1 constant omega:
eftcamb_params.update({'EFTCAMB_do_QSA': False})
pars_EFT = camb.set_params(
    **eftcamb_params, **bestfit, DoLateRadTruncation=True)

z = np.linspace(0., 2., 5)
# pars_EFT.set_matter_power(redshifts=z, kmax=1000.0);
# pars_EFT.NonLinear = model.NonLinear_none

pars_QSA.WantTransfer = True
pars_EFT.WantTransfer = True
###############################################################################
# run parameters:
results_QSA = camb.get_results(pars_QSA)
results_EFT = camb.get_results(pars_EFT)

k = np.array([0.01, 0.1, 1.0, 100, 1000])
#k= [1000.]
a = np.logspace(-10, 0, 2000)
vars_QSA, val_QSA = pars_QSA.EFTCAMB.get_scale_evolution(results_QSA, k, a)
vars_eft, val_eft = pars_EFT.EFTCAMB.get_scale_evolution(results_EFT, k, a)
#
out_root = out_dir+'test'


Hcrossind = np.argmin(np.abs((val_eft['adotoa']-k[:, np.newaxis]*np.ones(
    (len(k), 2000)))/k[:, np.newaxis]*np.ones((len(k), 2000))), axis=1)
Hcross = a[Hcrossind]


plt.figure(figsize=(5, 4))
for ind, _k in enumerate(k):

    plt.loglog(a, val_eft['QS_sigma_gen'][ind, :],
                              lw=1., ls='-', label="QSA, given full run")
    plt.loglog(a, val_QSA['QS_sigma_gen'][ind, :],
                              lw=1., ls='-', label="QSA, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("QSA $\Sigma$")
    plt.tight_layout()
    plt.show()


    plt.loglog(a, val_eft['QS_mu_gen'][ind, :],
                              lw=1., ls='-', label="QSA, given full run")
    plt.loglog(a, val_QSA['QS_mu_gen'][ind, :],
                              lw=1., ls='-', label="QSA, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("QSA $\mu$")
    plt.tight_layout()
    plt.show()

    plt.loglog(a, abs(val_eft['QS_pi'][ind, :]),
                              lw=1., ls='-', label="QSA, given full run")
    plt.loglog(a, abs(val_QSA['QS_pi'][ind, :]),
                              lw=1., ls='-', label="QSA, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("QSA $\pi$")
    plt.tight_layout()
    plt.show()


    plt.loglog(a, abs(val_eft['QS_pidot'][ind, :]),
                              lw=1., ls='-', label="QSA, given full run")
    plt.loglog(a, abs(val_QSA['QS_pidot'][ind, :]),
                              lw=1., ls='-', label="QSA, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("QSA $\dot{\pi}$")
    plt.tight_layout()
    plt.show()

plt.figure(figsize=(5, 4))
for ind, _k in enumerate(k):

    plt.loglog(a, val_eft['sigma_eft'][ind, :],
               lw=1., ls='-', label="full, given full run")
    plt.loglog(a, val_QSA['sigma_eft'][ind, :],
               lw=1., ls='-', label="full, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("$\Sigma$")
    plt.tight_layout()
    plt.show()


    plt.loglog(a, val_eft['mu'][ind, :],
               lw=1., ls='-', label="full, given full run")
    plt.loglog(a, val_QSA['mu'][ind, :],
               lw=1., ls='-', label="full, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("$\mu$")
    plt.tight_layout()
    plt.show()

    if (np.isnan(val_QSA['pi'][ind, :]).any()):
        plt.title("QSA run pi contains nans")
    plt.loglog(a, abs(val_eft['pi'][ind, :]),
               lw=1., ls='-', label="full, given full run")
    plt.loglog(a, abs(val_QSA['pi'][ind, :]),
               lw=1., ls='-', label="full, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("$\pi$")
    plt.tight_layout()
    plt.show()


    plt.loglog(a, abs(val_eft['pidot'][ind, :]),
               lw=1., ls='-', label="full, given full run")
    plt.loglog(a, abs(val_QSA['pidot'][ind, :]),
               lw=1., ls='-', label="full, given QSA run")
    plt.axvline(x=Hcross[ind], c='red', label="Hcross")
    plt.legend()
    plt.xlabel("$a$")
    plt.ylabel("$\dot{\pi}$")
    plt.tight_layout()
    plt.show()


print("trying to get matter power spectrum")
PK_QSA = camb.get_matter_power_interpolator(pars_QSA, hubble_units=False, k_hunit=False, kmax=10.0, zmax=10,nonlinear=True) #, extrap_kmax= 10**10) #,log_interp=False)
#PK_EFT = camb.get_matter_power_interpolator(pars_EFT, hubble_units=False, k_hunit=False, kmax=1000.0, zmax=10,nonlinear=True, extrap_kmax= 10**10) #,log_interp=False)
k = np.logspace(-5, 5,1000)
#plt.loglog(k, PK_EFT.P(0, k), label = 'full')
plt.loglog(k, PK_QSA.P(0, k), label = 'QSA')
plt.ylabel("P(k)")
plt.xlabel("k")
plt.legend()
plt.show()
