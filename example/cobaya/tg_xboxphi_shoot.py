from typing import Iterable, Sequence, Tuple, Union
from cobaya.typing import InfoDict
import numpy as np
from scipy.optimize import root, root_scalar
from scipy.special import lambertw

from cobaya.theory import Theory
from cobaya.log import LoggedError
from cobaya.yaml import yaml_load_file
from cobaya.run import run
from cobaya.theories.cosmo.boltzmannbase import PowerSpectrumInterpolator
from cobaya.conventions import Const
from cobaya.likelihood import Likelihood

import sys

class tg_xboxphi_shoot(Theory):
    cosmo_path = None
    extra_args = {}

    def initialize(self):
        if self.cosmo_path != None:
            sys.path.insert(1, self.cosmo_path)
        import camb
        from camb import CAMBError
        self.cosmo_shooter = camb
        self.cosmo_error = CAMBError
        
        self.extargs = {'dark_energy_model': 'EFTCAMB',
                        'EFTflag' :5,
                        'Horndeski_model': 7, 
                        'Horndeski_freefunc0_model': 11,
                        'Horndeski_freefunc1_model': 0,
                        'Horndeski_freefunc2_model': 1,
                        'Horndeski_freefunc3_model': 3,
                        'Horndeskic3Exp': 2,
                        'Horndeski_freefunc4_model': 0,
                        'Horndeski_freefunc5_model': 0,
                        'Horndeski_freefunc6_model': 0,
                        'Horndeski_parameter_number': 0,
                        'Horndeski_model_specific_ic': False,
                        'Horndeski_phidot_ini': 0,
                        'Horndeski_evolve_hubble': False,
                        'Horndeski_shooting': True,
                        'model_background_num_points': 10000,
                        'EFTCAMB_skip_stability': True,
                        'EFTCAMB_skip_RGR': True,
                        'model_background_a_ini': 1e-14,
                        'EFTCAMB_use_background': True,
                        'num_massive_neutrinos': 1,
                        'mnu': 0.06,
                        'nnu': 3.044,
                        }
        self.extargs.update(self.extra_args)
        self.cosmo_pars = {'ombh2': None, 'omch2': None, 'H0': None, 'Horndeskic01': None, 'Horndeskic20': None, 'Horndeskic30': None, 'Horndeski_phi_ini': None, }

    def get_can_provide_params(self):
        ps = ['Horndeskic00']
        return ps
    
    def get_requirements(self):
        rqs = self.cosmo_pars.copy()
        return rqs

    def calculate(self, state, want_derived=True, **params_values_dict):
        for par in self.cosmo_pars.keys():
            self.cosmo_pars[par] = self.provider.get_param(par)

        self.h0 = self.provider.get_param('H0')
        self.h0_mpc = 3.33564e-6*self.h0

        try:
            # rlt = root(self.is_ic, x0, tol=1e-3, options={'eps': 0.01})
            rlt = root_scalar(self.is_ic_scalar, rtol=1e-3, bracket=[0.01,1], x0=0.7, method="brentq")
        except self.cosmo_error:
            return False
        # if not rlt.success:
        #     self.log.debug("Shooting for horndeski parameters failed. "
        #                    "Assigning 0 likelihood and going on. ")
        #     return False
        # v = rlt.x[0]

        if not rlt.converged:
            self.log.debug("Shooting for horndeski parameters failed. "
                           "Assigning 0 likelihood and going on. ")
            return False
        v = rlt.root

        v0 = 3*self.h0_mpc**2*v

        state['derived'] = {'Horndeskic00': v0}

        return True

    def is_ic(self, x):
        v = x[0]
        v0 = 3*self.h0_mpc**2*v
        locparams = self.extargs.copy()
        locparams.update({'Horndeskic00': v0})
        locparams.update(self.cosmo_pars)
        pars = self.cosmo_shooter.set_params(**locparams)
        res = self.cosmo_shooter.get_background(pars, no_thermo=True)
        rlt = np.empty(1)
        rlt[0] = res.Params.EFTCAMB_parameter_cache.h0/self.h0 - 1
        return rlt
    
    def is_ic_scalar(self, x):
        v = x
        v0 = 3*self.h0_mpc**2*v
        locparams = self.extargs.copy()
        locparams.update({'Horndeskic00': v0})
        locparams.update(self.cosmo_pars)
        pars = self.cosmo_shooter.set_params(**locparams)
        res = self.cosmo_shooter.get_background(pars, no_thermo=True)
        rlt = res.Params.EFTCAMB_parameter_cache.h0/self.h0 - 1
        return rlt

