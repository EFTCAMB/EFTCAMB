"""
This file contains the basic tools to interface EFTCAMB to camb.

@author : Fabrizio Renzi, Marco Raveri

"""

import re
import warnings
from ctypes import POINTER, addressof, byref, c_bool, c_double, c_int, c_size_t, create_string_buffer

import numpy as np

from ._config import config
from .baseconfig import AllocatableArrayDouble, CAMBValueError, F2003Class, camblib, fortran_class

###############################################################################
# Interface functions:
###############################################################################


def _import_string_function(bind_c_name, module_name):
    """
    Import a Fortran routine with character arguments.

    The module procedures take character(len=200) dummies, for which every supported compiler (gfortran, ifort,
    flang) appends a hidden length argument (size_t, by value) after the explicit arguments, in argument order.
    A bind(c) version exported under bind_c_name, with the same explicit lengths in the same positions, is used
    instead when the library provides one.
    """
    try:
        return camblib[bind_c_name]
    except AttributeError:
        return getattr(camblib, module_name)


# length of the character dummies of the Fortran getters and setters
_STRING_LENGTH = 200
# explicit value of the hidden length argument of those character dummies
_STRING_LENGTH_ARG = c_size_t(_STRING_LENGTH)


def _string_buffer(value=""):
    """
    Writable, blank padded buffer for a character(len=200) argument.
    Longer values are truncated to 200 bytes on a character boundary.
    """
    data = value.encode("utf-8")[:_STRING_LENGTH].decode("utf-8", errors="ignore").encode("utf-8")
    return create_string_buffer(data.ljust(_STRING_LENGTH), _STRING_LENGTH)


def _string_value(buffer):
    """
    Python string held by a buffer filled by Fortran.
    """
    return buffer.raw.decode("utf-8", errors="replace").replace("\x00", " ").strip()


# Fortran real literal with a d/D exponent (e.g. 1.d-4), which Python float() does not accept
_FORTRAN_D_REAL = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)[dD][+-]?\d+")


# ini keys of the EFTCAMB flags read at initialization, where they differ from the python field names
_FLAG_INI_KEYS = {"EFTCAMB_feedback_level": "feedback_level", "EFTCAMB_pert_turn_on": "EFTCAMB_turn_on_time"}


def _format_ini_value(value):
    """
    Ini string of a python bool, int, float or string value.
    """
    if isinstance(value, (bool, np.bool_)):
        return "T" if value else "F"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _convert_ini_value(value):
    """
    Convert a string read from an ini file to bool ('T'/'F'), int or float when possible, otherwise keep the string.
    """
    if value == "T":
        return True
    if value == "F":
        return False
    for fn in (int, float):
        try:
            return fn(value)
        except ValueError:
            pass
    if _FORTRAN_D_REAL.fullmatch(value):
        return float(value.replace("d", "e").replace("D", "e"))
    return value


# Names of EFTCAMB parameters, used to accept (and ignore) them in camb.set_params when EFTCAMB is off (EFTflag=0).
# The keys of the model parameters are partly built at run time (parametrization names), so they are recognised by
# their prefixes; this is not an exhaustive list.
_EFTCAMB_KEY_PREFIXES = (
    "EFT",
    "PureEFT",
    "RPH",
    "OL",
    "Designer3p1",
    "Horava",
    "Horndeski",
    "Beyond_Horndeski",
    "ExtendedGalileon",
    "Scaling_Cubic",
    "Shift_Symmetric",
    "Kmimic",
)
_EFTCAMB_KEYS = frozenset(
    (
        "AltParEFTmodel",
        "DesignerEFTmodel",
        "FullMappingEFTmodel",
        "feedback_level",
        "output_root",
        "alphaU",
        "gammaA",
        "gammaU",
        "drag_initial_conditions",
        "eps2_0",
        "f_ac",
        "Log_ac",
        "field_max",
        "field_min",
        "m",
        "p",
        "wf",
        "phidot_ini",
        "potential_model",
        "model_background_a_ini",
        "model_background_a_final",
        "model_background_num_points",
        "want_debug",
        "want_LCDM",
    )
)


def is_eftcamb_parameter(key):
    """
    Whether key looks like the name of an EFTCAMB parameter (flag or model parameter).
    """
    return key in _EFTCAMB_KEYS or key.startswith(_EFTCAMB_KEY_PREFIXES)


EFTCAMB_initialize = camblib.__eftcamb_python_MOD_initialize_eftcamb_from_py
# bind(c) routine: its name is fixed by the binding label for every compiler
EFTCAMB_initialize_read_params = camblib.initialize_EFTCAMB_from_py_read_params_

EFTCAMB_feedback = camblib.__eftcamb_python_MOD_eftcamb_feedback
EFTCAMB_get_read_parameters = _import_string_function(
    "eftcamb_get_read_parameters_c", "__eftcamb_python_MOD_get_read_parameters"
)

# EFTCAMB getters:
EFTCAMB_get_num_params = camblib.__eftcamb_python_MOD_eftcamb_get_num_params
EFTCAMB_get_param_names = _import_string_function(
    "eftcamb_get_param_names_c", "__eftcamb_python_MOD_eftcamb_get_param_names"
)
EFTCAMB_get_param_labels = _import_string_function(
    "eftcamb_get_param_labels_c", "__eftcamb_python_MOD_eftcamb_get_param_labels"
)
EFTCAMB_get_param_values = camblib.__eftcamb_python_MOD_eftcamb_get_param_values
EFTCAMB_get_model_name = _import_string_function(
    "eftcamb_get_model_name_c", "__eftcamb_python_MOD_eftcamb_get_model_name"
)
EFTCAMB_get_effective_w0wa_values = camblib.__eftcamb_python_MOD_eftcamb_get_effective_w0wa_values

# EFTCAMB setters:
EFTCAMB_set_model_name = _import_string_function(
    "eftcamb_set_model_name_c", "__eftcamb_python_MOD_eftcamb_set_model_name"
)

# EFTCAMB printing variables:
EFTCAMB_quantities_evolution = camblib.__eftcamb_python_MOD_eftcamb_getevolution

# EFTCAMB printing EFT functions (no recomputation):
EFTCAMB_eftfunctions = camblib.__eftcamb_python_MOD_eftcamb_geteftfunctions

###############################################################################
# EFTCAMB classes:
###############################################################################


@fortran_class
class EFTCAMB(F2003Class):
    _fortran_class_module_ = "EFTCAMB_main"
    _fortran_class_name_ = "TEFTCAMB"

    _fields_ = (
        # Model selection flags:
        ("EFTflag", c_int, "Main EFTCAMB model selection flag. Decides one of the four modes to run EFTCAMB."),
        ("PureEFTmodel", c_int, "Model selection flag for pure EFT models"),
        ("AltParEFTmodel", c_int, "Model selection flag for alternative EFT parameterization"),
        ("DesignerEFTmodel", c_int, "Model selection flag for designer EFT models"),
        ("FullMappingEFTmodel", c_int, "Model selection flag for full mapping EFT models"),
        # Stability flags
        (
            "EFT_ghost_math_stability",
            c_bool,
            "Flag that decides whether to enforce the mathematical ghost condition. This is less conservative than the physical requirement. Enforces the condition on both scalar and tensor sectors",
        ),
        (
            "EFT_mass_math_stability",
            c_bool,
            "Flag that decides whether to enforce the mathematical mass condition. This should prevent fast growing instabilities but is not fully rigorous. Use is deprecated",
        ),
        (
            "EFT_ghost_stability",
            c_bool,
            "Flag that decides whether to enforce the physical ghost condition. Works up to Horndeski.",
        ),
        (
            "EFT_gradient_stability",
            c_bool,
            "Flag that decides whether to enforce the physical gradient condition. Works up to Horndeski.",
        ),
        (
            "EFT_mass_stability",
            c_bool,
            "Flag that decides whether to enforce the physical mass condition. Works up to Horndeski.",
        ),
        (
            "EFT_mass_stability_rate",
            c_double,
            "Flag that sets the rate for the mass instability in units of Hubble time.",
        ),
        (
            "EFT_additional_priors",
            c_bool,
            "Flag that extablishes whether to use additional priors that are related to the specific model. Each model has the possibility of implementing their own additional stability requirements.",
        ),
        # Positvity bounds flags
        (
            "EFT_positivity_bounds",
            c_bool,
            "Flag that decides whether to enforce the positivity bounds. Works up to Horndeski. ",
        ),
        (
            "EFT_minkowski_limit",
            c_bool,
            "Flag that checks the existence of a healthy Minkowski limit. Works up to Horndeski. ",
        ),
        # QSA flags
        ("EFT_check_QSA", c_bool, "Flag that decides whether check QSA of mu"),
        ("EFT_QSA_threshold", c_double, "Tolerence for relative difference between (mu,sigma) and that of QSA"),
        ("EFT_QSA_time", c_double, "Minimum scale factor at which the code checks for QSA"),
        ("EFT_QSA_k", c_double, "Minimum k at which the code checks for QSA "),
        # Initial condition type
        (
            "EFT_IC_type",
            c_int,
            "Initial condition type. 1: GR, 2: assume Omega=const. at IC time, 3: assume Omega=const. Gamma=const. at IC time",
        ),
        # Working flags
        ("EFTCAMB_feedback_level", c_int, "Amount of feedback that is printed to screen"),
        (
            "EFTCAMB_back_turn_on",
            c_double,
            "Smallest scale factor that the code should consider when computing the background. Set to zero (default), change background at all times.",
        ),
        (
            "EFTCAMB_pert_turn_on",
            c_double,
            "Smallest scale factor at which EFTCAMB evolves perturbations. Default set to a=0.01.",
        ),
        (
            "EFTCAMB_GR_threshold",
            c_double,
            "Tollerance on the dimensionless EFT functions to be considered GR. Default is 10^-8.",
        ),
        ("EFTCAMB_stability_time", c_double, "Minimum scale factor at which the code checks for stability of a theory"),
        ("EFTCAMB_stability_threshold", c_double, "Threshold for the stability module to consider the model stable."),
        ("EFTCAMB_model_is_designer", c_bool, "Logical flag that establishes whether the model is designer or not."),
        (
            "EFTCAMB_effective_w0wa",
            c_bool,
            "Logical flag that establishes whether the model has effective w0wa or not.",
        ),
        (
            "EFTCAMB_use_background",
            c_bool,
            "Whether use the background conformal time, cosmic time and sound horizon computed by EFTCAMB instead of integrating them in CAMB. Currently only supported by the Horndeski module.",
        ),
        (
            "EFTCAMB_evolve_delta_phi",
            c_bool,
            "Flag that decides whether integrate delta_phi instead of pi for covariant theories",
        ),
        ("EFTCAMB_evolve_metric_h", c_bool, "Flag that decides whether integrate h' instead of using constraints"),
        ("EFTCAMB_skip_stability", c_bool, "Flag that decides whether skip all stability checks."),
        ("EFTCAMB_skip_RGR", c_bool, "Flag that decides whether skip the return to GR checks."),
        (
            "EFTCAMB_write_background",
            c_bool,
            "Flag that decides whether to write the background quantities when called from an ini file.",
        ),
    )

    def initialize_parameters(self, camb_parameters, EFTCAMB_params, print_header=True):
        """
        Function to initialize EFTCAMB given a dictionary with all the options.

        :param camb_parameters: the input CAMBparams class containing EFTCAMB.
        :param EFTCAMB_params: a dictionary containing all the EFTCAMB options.
        """
        # convert dictionary to list of strings that EFTCAMB can digest:
        _keys, _values = [], []
        for _key in sorted(EFTCAMB_params.keys()):
            # append key:
            _keys.append(_key)
            # append value:
            _val = EFTCAMB_params[_key]
            if isinstance(_val, (bool, np.bool_)):
                if _val:
                    _values.append("T")
                else:
                    _values.append("F")
            elif isinstance(_val, bytes):
                _values.append(_val.decode("utf-8"))
            else:
                _values.append(str(_val))
        # format the strings:
        _keys = [_key.encode("utf-8") for _key in _keys]
        _values = [_val.encode("utf-8") for _val in _values]
        # first call to initializer (this initializes the internal file):
        EFTCAMB_initialize_read_params(
            create_string_buffer(b" "), c_int(1), create_string_buffer(b" "), c_int(1), c_int(0)
        )
        # append to file (this adds the parameters to the file):
        for _key, _val in zip(_keys, _values):
            EFTCAMB_initialize_read_params(
                create_string_buffer(_key), c_int(len(_key)), create_string_buffer(_val), c_int(len(_val)), c_int(1)
            )
        # finish initialization (this calls the allocation of models):
        temp = c_int(0)
        EFTCAMB_initialize(byref(camb_parameters), byref(temp))
        if temp.value > 0:
            raise CAMBValueError("EFTCAMB error during initialization (see the printed output for the reason)")

        # the EFTCAMB object of camb_parameters is re-initialized in place, so handles to it stay valid.
        # A handle to the EFTCAMB object of other parameters does not see this initialization:
        eftcamb = camb_parameters.EFTCAMB
        if addressof(eftcamb) != addressof(self):
            warnings.warn(
                "initialize_parameters was called on the EFTCAMB object of other parameters, "
                "which is not changed. Use camb_parameters.EFTCAMB to access the initialized one.",
                stacklevel=2,
            )
        # print feedback:
        if print_header:
            eftcamb.feedback(print_params=True)

    def read_parameters(self, raw=False):
        """
        Function that returns all the parameters that have been read to
        initialize this EFTCAMB object.

        :param raw: if True, return the values as the strings that were read (as in an ini file), otherwise convert
            them to bool ('T'/'F'), int or float (including Fortran 'd' exponents, e.g. '1.d-4') when possible.
        :return: dictionary of parameter names and values, which can be passed back to
            :meth:`initialize_parameters` (or written to an ini file when raw=True)
        """
        # get the total number of parameters
        _key, _value = _string_buffer(), _string_buffer()
        num_params = c_int(0)
        EFTCAMB_get_read_parameters(
            byref(self), _key, _value, byref(num_params), byref(c_int(0)), _STRING_LENGTH_ARG, _STRING_LENGTH_ARG
        )
        num_params = num_params.value
        # get the parameters:
        params = {}
        for ind in range(num_params):
            index = c_int(ind + 1)
            _key, _value = _string_buffer(), _string_buffer()
            EFTCAMB_get_read_parameters(
                byref(self), _key, _value, byref(index), byref(c_int(1)), _STRING_LENGTH_ARG, _STRING_LENGTH_ARG
            )
            params[_string_value(_key)] = _string_value(_value)
        if raw:
            return params
        return {key: _convert_ini_value(value) for key, value in params.items()}

    def current_parameters(self):
        """
        Parameters (as ini strings) that re-create this EFTCAMB object with :meth:`initialize_parameters`: the
        parameters read at initialization, with the flags changed from python since then set to their current values.

        :return: dictionary of parameter names and string values
        """
        if not self.EFTflag:
            return {"EFTflag": "0"}
        params = self.read_parameters(raw=True)
        for field, _ in self.get_all_fields():
            key = _FLAG_INI_KEYS.get(field, field)
            if key in params:
                value = getattr(self, field)
                if _convert_ini_value(params[key]) != value:
                    params[key] = _format_ini_value(value)
        return params

    def feedback(self, print_params=True):
        """
        Feedback for EFTCAMB.

        :param print_params: (optional) logical parameter that decides whether
            to print the values of model parameters or not.
        """
        # logical (default kind) argument: passed as a C int, like CAMB's logical fields
        EFTCAMB_feedback(byref(self), byref(c_int(1 if print_params else 0)))

    def num_params(self):
        """
        Returns the number of additional parameters of the EFTCAMB model.
        """
        temp = c_int(0)
        EFTCAMB_get_num_params(byref(self), byref(temp))
        return temp.value

    def param_names(self):
        """
        Return the parameter names for the EFTCAMB model.
        """
        n_params = self.num_params()
        param_names = []
        for i in range(n_params):
            name = _string_buffer()
            EFTCAMB_get_param_names(byref(self), byref(c_int(i + 1)), name, _STRING_LENGTH_ARG)
            param_names.append(_string_value(name))
        return param_names

    def param_labels(self):
        """
        Return the parameter labels for the EFTCAMB model.
        """
        # get number of parameters:
        n_params = self.num_params()
        # get code labels:
        param_labels = []
        for i in range(n_params):
            name = _string_buffer()
            EFTCAMB_get_param_labels(byref(self), byref(c_int(i + 1)), name, _STRING_LENGTH_ARG)
            param_labels.append(_string_value(name))
        return param_labels

    def param_values(self):
        """
        Return the parameter values for the EFTCAMB model being used.
        """
        n_params = self.num_params()
        param_values = []
        for i in range(n_params):
            temp = c_double(0)
            EFTCAMB_get_param_values(byref(self), byref(c_int(i + 1)), byref(temp))
            param_values.append(temp.value)
        return param_values

    def model_name(self):
        """
        Return the model name for the EFTCAMB model.
        """
        model = _string_buffer()
        EFTCAMB_get_model_name(byref(self), model, _STRING_LENGTH_ARG)
        return _string_value(model)

    def set_model_name(self, name):
        """
        Set the model name for the EFTCAMB model.

        :param name: string with the new name of the model (at most 200 bytes in UTF-8 are used).
        """
        EFTCAMB_set_model_name(byref(self), _string_buffer(name), _STRING_LENGTH_ARG)

    def get_effective_w0wa(self):
        """
        Returns the effective w0 wa values.
        """
        w0 = c_double(0)
        wa = c_double(0)
        EFTCAMB_get_effective_w0wa_values(byref(self), byref(w0), byref(wa))
        return w0.value, wa.value

    def get_scale_evolution(self, camb_results, q, a):
        """
        Get the time evolution of the input EFTCAMB model at the requested scale factors.

        :param camb_results: a :class:`~.results.CAMBdata`
        :param q: wavenumber values to calculate (or array of k values)
        :param a: array of requested scale factors to output (any order)
        :return: as :meth:`get_evolution`, with the times in the order of the input scale factors
        """
        # convert scale factor to redshift:
        z = 1.0 / np.atleast_1d(np.asarray(a, dtype=np.float64)) - 1.0
        # then to conformal time:
        eta = camb_results.conformal_time(z)
        # then call the other function:
        return self.get_evolution(camb_results, q, eta)

    def get_evolution(self, camb_results, q, eta):
        """
        Get the time evolution of the input EFTCAMB model.

        This function is redoing calculations so it should be used only for pretty plots.

        :param camb_results: a :class:`~.results.CAMBdata`
        :param q: wavenumber values to calculate (or array of k values)
        :param eta: array of requested conformal times to output (any order)
        :return: list of field names, and structured array of timestep cache values with shape (len(q), len(eta)),
            in the order of the input wavenumbers and times
        """
        # digest the input scale array:
        if np.isscalar(q):
            k = np.array([q], dtype=np.float64)
        else:
            k = np.array(q, dtype=np.float64)
        k = np.asfortranarray(k)
        # digest the input time array:
        times = np.asarray(np.atleast_1d(eta), dtype=np.float64)
        # the Fortran integration needs increasing times: sort them and restore the caller's order on output
        indices = np.argsort(times)
        sorted_times = np.asfortranarray(times[indices])
        # prepare the output array with the time step cache type:
        results = (EFTCAMB_timestep_cache * times.shape[0] * k.shape[0])()
        # call the fortran function to do the heavy lifting:
        EFTCAMB_quantities_evolution(
            byref(camb_results),
            byref(c_int(k.shape[0])),
            k.ctypes.data_as(POINTER(c_double)),
            byref(c_int(times.shape[0])),
            sorted_times.ctypes.data_as(POINTER(c_double)),
            byref(results),
        )
        config.check_global_error("get_evolution")
        # now unpackage the results, as [k][time] in the order of the input times:
        fields = results[0][0].variables()
        results = np.ctypeslib.as_array(results)
        i_rev = np.zeros(times.shape, dtype=int)
        i_rev[indices] = np.arange(times.shape[0])
        return fields, results[:, i_rev]

    def get_eft_functions(self, camb_results, a):
        """
        Get the EFT functions (those depends only on the background) of the input EFTCAMB model.

        :param camb_results: a :class:`~.results.CAMBdata`
        :param a: array of requested scale factors to output
        """
        # convert scale factor to redshift:
        if np.isscalar(a):
            z = 1.0 / np.array([a], dtype=np.float64) - 1.0
        else:
            z = 1.0 / np.array(a, dtype=np.float64) - 1.0
        # then to conformal time:
        eta = camb_results.conformal_time(z)
        # digest the input scale factor array:
        times = np.asfortranarray(np.atleast_1d(a), dtype=np.float64)
        etas = np.asfortranarray(np.atleast_1d(eta), dtype=np.float64)
        # prepare the output array with the time step cache type:
        results = (EFTCAMB_timestep_cache * times.shape[0])()
        # call the fortran function:
        config.global_error_flag = 0
        EFTCAMB_eftfunctions(
            byref(camb_results),
            byref(c_int(times.shape[0])),
            times.ctypes.data_as(POINTER(c_double)),
            etas.ctypes.data_as(POINTER(c_double)),
            byref(results),
        )
        config.check_global_error("get_eft_functions")
        fields = results[0].variables()
        results = np.ctypeslib.as_array(results)
        return fields, results


@fortran_class
class EFTCAMB_parameter_cache(F2003Class):
    _fortran_class_module_ = "EFTCAMB_cache"
    _fortran_class_name_ = "TEFTCAMB_parameter_cache"

    _fields_ = (
        # 1) relative densities:
        ("omegac", c_double),
        ("omegab", c_double),
        ("omegav", c_double),
        ("omegak", c_double),
        ("omegan", c_double),
        ("omegag", c_double),
        ("omegar", c_double),
        # 2) Hubble constant:
        ("h0", c_double),
        ("h0_Mpc", c_double),
        # 3) densities:
        ("grhog", c_double),
        ("grhornomass", c_double),
        ("grhoc", c_double),
        ("grhob", c_double),
        ("grhov", c_double),
        ("grhok", c_double),
        # 4) massive neutrinos:
        ("Num_Nu_Massive", c_int),
        ("Nu_mass_eigenstates", c_int),
        ("grhormass", AllocatableArrayDouble),
        ("nu_masses", AllocatableArrayDouble),
        # 5a) Horndeski module: EDE shooting parameters
        ("maxfrac_hdsk", c_double),
        ("z_maxfrac_hdsk", c_double),
        # 5b) Horndeski module: designer approach products
    )


@fortran_class
class EFTCAMB_timestep_cache(F2003Class):
    """
    This class contains the time step cache that is internally used by EFTCAMB for all calculations.
    Right now the fields must be updated if the fortran code is modified.
    In the future we can build this interface automatically.
    """

    _fortran_class_module_ = "EFTCAMB_cache"
    _fortran_class_name_ = "TEFTCAMB_timestep_cache"

    _fields_ = (
        # 1) time and k:
        ("a", c_double),  # the value of the scale factor at which the cache is being used.
        ("tau", c_double),  # the value of conformal time at the given scale factor.
        ("k", c_double),  # the scale that is being solved for. In \f$ \mbox{Mpc}^{-1} \f$.
        # 2) total matter densities:
        ("grhoa2", c_double),  # the input value of \f$ \sum_m\rho_m / a^2 m_0^2 \f$.
        ("grhom_t", c_double),  # the value of \f$ \sum_m\rho_m a^2 /m_0^2 \f$.
        ("gpresm_t", c_double),  # the value of \f$ \sum_m P_m a^2 /m_0^2 \f$.
        ("gpresdotm_t", c_double),  # the value of \f$ \sum_m\dot{P}_m a^2 /m_0^2 \f$.
        ("gpresdotdotm_t", c_double),  # < the value of \f$ \sum_m \ddot{P}_m a^2 /m_0^2 \f$.
        # 3) densities and pressure of the various species:
        ("grhob_t", c_double),  # the value of \f$ \rho_b a^2 / m_0^2 \f$
        ("grhoc_t", c_double),  # the value of \f$ \rho_{cdm} a^2 / m_0^2 \f$
        ("grhor_t", c_double),  # the value of \f$ \rho_{\nu} a^2 / m_0^2 \f$
        ("grhog_t", c_double),  # the value of \f$ \rho_{\gamma} a^2 / m_0^2 \f$
        (
            "grhov_t",
            c_double,
        ),  # the value of \f$ \rho_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        (
            "gpiv_t",
            c_double,
        ),  # the value of \f$ \P_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        ("grhonu_tot", c_double),  # the value of \f$ \sum_\nu \rho_{m\nu} a^2 / m_0^2 \f$
        ("gpinu_tot", c_double),  # the value of \f$ \sum_\nu P_{m\nu} a^2 / m_0^2 \f$
        ("grhonudot_tot", c_double),  # the value of \f$ \sum_\nu \dot{\rho}_{m\nu} a^2 / m_0^2 \f$
        ("gpinudot_tot", c_double),  # the value of \f$ \sum_\nu \dot{P}_{m\nu} a^2 / m_0^2 \f$
        ("gpinudotdot_tot", c_double),  # < the value of \f$ \sum_\nu \ddot{P}_{m\nu} a^2 / m_0^2 \f$
        # 4) expansion history:
        ("adotoa", c_double),  # the value of \f$ \mathcal{H} \f$ at the given scale factor.
        ("Hdot", c_double),  # the value of \f$ d\mathcal{H} /d \tau \f$ at the given scale factor.
        ("Hdotdot", c_double),  # the value of \f$ d^2 \mathcal{H} / d \tau^2 \f$ at the given scale factor.
        ("Hdotdotdot", c_double),  # < the value of \f$ d^3 \mathcal{H} / d \tau^3 \f$ at the given scale factor.
        # 5) EFT functions:
        ("EFTOmegaV", c_double),  # the value of Omega \f$ \Omega(a) \f$.
        ("EFTOmegaP", c_double),  # the value of the derivative wrt scale factor of Omega \f$ d \Omega(a) / da \f$.
        (
            "EFTOmegaPP",
            c_double,
        ),  # the value of the second derivative wrt scale factor of Omega \f$ d^2 \Omega(a) / da^2 \f$.
        (
            "EFTOmegaPPP",
            c_double,
        ),  # the value of the third derivative wrt scale factor of Omega \f$ d^3 \Omega(a) / da^3 \f$.
        (
            "EFTOmegaPPPP",
            c_double,
        ),  # < the value of the fourth derivative wrt scale factor of Omega \f$ d^4 \Omega(a) / da^4 \f$.
        ("EFTc", c_double),  # the value of \f$ c a^2/m_0^2 \f$.
        ("EFTcdot", c_double),  # the value of \f$ \dot{c} a^2/m_0^2 \f$.
        ("EFTcdotdot", c_double),  # the value of \f$ \dot{c} a^2/m_0^2 \f$.
        ("EFTLambda", c_double),  # the value of \f$ \Lambda a^2/m_0^2 \f$.
        (
            "EFTLambdadot",
            c_double,
        ),  # the value of \f$ \dot{\Lambda}a^2/m_0^2 \f$. Derivative of \f$ \Lambda\f$ wrt conformal time.
        (
            "EFTLambdadotdot",
            c_double,
        ),  # < the value of \f$ \ddot{\Lambda}a^2/m_0^2 \f$. Derivative of \f$ \Lambda\f$ wrt conformal time.
        ("EFTGamma1V", c_double),  # the value of Gamma 1 \f$ \gamma_1(a) \f$.
        (
            "EFTGamma1P",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 1 \f$  d \gamma_1(a) / da \f$.
        (
            "EFTGamma1PP",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 1 \f$  d \gamma_1(a) / da \f$.
        ("EFTGamma2V", c_double),  # the value of Gamma 2 \f$ \gamma_2(a) \f$.
        (
            "EFTGamma2P",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 2 \f$  d \gamma_2(a) / da \f$.
        (
            "EFTGamma2PP",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 2 \f$  d \gamma_2(a) / da \f$.
        (
            "EFTGamma2PPP",
            c_double,
        ),  # < the value of the third derivative wrt scale factor of Gamma 3 \f$ d^3 \gamma_2(a) / da^3 \f$.
        ("EFTGamma3V", c_double),  # the value of Gamma 3 \f$ \gamma_3(a) \f$.
        (
            "EFTGamma3P",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 3 \f$  d \gamma_3(a) / da \f$.
        (
            "EFTGamma3PP",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 3 \f$  d \gamma_3(a) / da \f$.
        (
            "EFTGamma3PPP",
            c_double,
        ),  # < the value of the third derivative wrt scale factor of Gamma 3 \f$ d^3 \gamma_3(a) / da^3 \f$.
        (
            "EFTGamma3PPPP",
            c_double,
        ),  # < the value of the fourth derivative wrt scale factor of Gamma 3 \f$ d^4 \gamma_3(a) / da^4 \f$.
        ("EFTGamma4V", c_double),  # the value of Gamma 4 \f$ \gamma_4(a) \f$.
        (
            "EFTGamma4P",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 4 \f$  d \gamma_4(a) / da \f$.
        (
            "EFTGamma4PP",
            c_double,
        ),  # the value of the second derivative wrt scale factor of Gamma 4 \f$  d^2 \gamma_4(a) / da^2 \f$.
        ("EFTGamma5V", c_double),  # the value of Gamma 5 \f$ \gamma_5(a) \f$.
        (
            "EFTGamma5P",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 5 \f$  d \gamma_5(a) / da \f$.
        ("EFTGamma6V", c_double),  # the value of Gamma 6 \f$ \gamma_6(a) \f$.
        (
            "EFTGamma6P",
            c_double,
        ),  # the value of the derivative wrt scale factor of Gamma 6 \f$  d \gamma_6(a) / da \f$.
        # 5b) alternative EFT functions
        ("alphaB", c_double),  # the value of the Horndeski function \f$ \alpha_B \f$
        ("alphaBdot", c_double),  # the value of the derivative wrt conformal time of alphaB \f$ \dot{\alpha}_B \f$
        ("alphaT", c_double),  # the value of the Horndeski function \f$ \alpha_T \f$
        ("alphaTdot", c_double),  # the value of the derivative wrt conformal time of alphaT \f$ \dot{\alpha}_T \f$
        ("alphaK", c_double),  # the value of the Horndeski function \f$ \alpha_K \f$
        ("alphaKdot", c_double),  # the value of the derivative wrt conformal time of alphaK \f$ \dot{\alpha}_K \f$
        ("alphaM", c_double),  # the value of the Horndeski function \f$ \alpha_M \f$
        ("alphaMdot", c_double),  # the value of the derivative wrt conformal time of alphaK \f$ \dot{\alpha}_M \f$
        ("Meff2", c_double),  # the value of the effective Planck mass over today's measured value \f$ M_*^2/m_0^2 \f$
        # 5c) delta phi perturbation coefficients
        (
            "spi1",
            c_double,
        ),  # the coefficient before delta_phi_prime on the LHS of the delta_phi equation, see the appendix of XXXX.XXXX
        (
            "spi2",
            c_double,
        ),  # the coefficient before delta_phi on the LHS of the delta_phi equation, see the appendix of XXXX.XXXX
        (
            "spi3",
            c_double,
        ),  # the coefficient before k^2*delta_phi on the LHS of the delta_phi equation, see the appendix of XXXX.XXXX
        (
            "spih",
            c_double,
        ),  # the coefficient before h_prime on the RHS of the delta_phi equation, see the appendix of XXXX.XXXX
        (
            "spie",
            c_double,
        ),  # the coefficient before k^2*eta on the RHS of the delta_phi equation, see the appendix of XXXX.XXXX
        (
            "spim",
            c_double,
        ),  # the coefficient before P_matter on the RHS of the delta_phi equation, see the appendix of XXXX.XXXX
        (
            "s00",
            c_double,
        ),  # the coefficient before delta_phi in the modified 00 einstein equation, see the appendix of XXXX.XXXX
        (
            "s00k",
            c_double,
        ),  # the coefficient before k^2*delta_phi in the modified 00 einstein equation, see the appendix of XXXX.XXXX
        (
            "s00p",
            c_double,
        ),  # the coefficient before delta_phi_prime in the modified 00 einstein equation, see the appendix of XXXX.XXXX
        (
            "s0i",
            c_double,
        ),  # the coefficient before delta_phi in the modified 0i einstein equation, see the appendix of XXXX.XXXX
        (
            "s0ip",
            c_double,
        ),  # the coefficient before delta_phi_prime in the modified 0i einstein equation, see the appendix of XXXX.XXXX
        (
            "sii",
            c_double,
        ),  # the coefficient before delta_phi in the modified ii einstein equation, see the appendix of XXXX.XXXX
        (
            "siik",
            c_double,
        ),  # the coefficient before k^2*delta_phi in the modified ii einstein equation, see the appendix of XXXX.XXXX
        (
            "siip",
            c_double,
        ),  # the coefficient before delta_phi_prime in the modified ii einstein equation, see the appendix of XXXX.XXXX
        (
            "siipp",
            c_double,
        ),  # the coefficient before delta_phi_prime_prime in the modified ii einstein equation, see the appendix of XXXX.XXXX
        (
            "sij",
            c_double,
        ),  # the coefficient before delta_phi in the modified ij einstein equation, see the appendix of XXXX.XXXX
        (
            "sijdot",
            c_double,
        ),  # the derivative wrt conformal time of the coefficient before delta_phi in the modified ij einstein equation, see the appendix of XXXX.XXXX
        # 5d) background scalar field
        ("phi_scf", c_double),  # the background value of the scalar field phi
        ("V_scf", c_double),  # the value of the scalar field potential
        ("dphi", c_double),  # the derivative of the scalar field phi wrt conformal time
        ("ddphi", c_double),  # the second derivative of the scalar field phi wrt conformal time
        ("dddphi", c_double),  # the third derivative of the scalar field phi wrt conformal time
        # 5e) scalar field time scale
        (
            "dtau_hdsk",
            c_double,
        ),  # the conformal time difference between two nearby points of phi=0, estimation of the field evolution timescale if it is oscillating
        # 6) other background quantities:
        (
            "grhoq",
            c_double,
        ),  # the value of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        (
            "gpresq",
            c_double,
        ),  # the value of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        (
            "grhodotq",
            c_double,
        ),  # the value of the time derivative of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        (
            "gpresdotq",
            c_double,
        ),  # the value of the time derivative of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        # 7) the Einstein equations coefficients:
        (
            "EFTeomF",
            c_double,
        ),  # the value of the Einstein equations coefficient F. Refer to the Numerical Notes for the definition.
        (
            "EFTeomN",
            c_double,
        ),  # the value of the Einstein equations coefficient N. Refer to the Numerical Notes for the definition.
        (
            "EFTeomNdot",
            c_double,
        ),  # the value of the Einstein equations coefficient dN/dtau. Refer to the Numerical Notes for the definition.
        (
            "EFTeomX",
            c_double,
        ),  # the value of the Einstein equations coefficient X. Refer to the Numerical Notes for the definition.
        (
            "EFTeomXdot",
            c_double,
        ),  # the value of the Einstein equations coefficient dX/dtau. Refer to the Numerical Notes for the definition.
        (
            "EFTeomY",
            c_double,
        ),  # the value of the Einstein equations coefficient Y. Refer to the Numerical Notes for the definition.
        (
            "EFTeomG",
            c_double,
        ),  # the value of the Einstein equations coefficient G. Refer to the Numerical Notes for the definition.
        (
            "EFTeomU",
            c_double,
        ),  # the value of the Einstein equations coefficient U. Refer to the Numerical Notes for the definition.
        (
            "EFTeomL",
            c_double,
        ),  # the value of the Einstein equations coefficient L. Refer to the Numerical Notes for the definition.
        (
            "EFTeomM",
            c_double,
        ),  # the value of the Einstein equations coefficient M. Refer to the Numerical Notes for the definition.
        (
            "EFTeomV",
            c_double,
        ),  # the value of the Einstein equations coefficient V. Refer to the Numerical Notes for the definition.
        (
            "EFTeomVdot",
            c_double,
        ),  # the value of the Einstein equations coefficient dV/dtau. Refer to the Numerical Notes for the definition.
        (
            "EFTeomQ",
            c_double,
        ),  # the value of the Einstein equations coefficient Q. Refer to the Numerical Notes for the definition.
        # 8) pi field factors:
        (
            "EFTpiA1",
            c_double,
        ),  # the value of the pi field equation coefficient A1. Scale independent part of A. Refer to the Numerical Notes for the definition.
        (
            "EFTpiA2",
            c_double,
        ),  # the value of the pi field equation coefficient A2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        (
            "EFTpiB1",
            c_double,
        ),  # the value of the pi field equation coefficient B1. Scale independent part of B. Refer to the Numerical Notes for the definition.
        (
            "EFTpiB2",
            c_double,
        ),  # the value of the pi field equation coefficient B2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        (
            "EFTpiC",
            c_double,
        ),  # the value of the pi field equation coefficient C. Refer to the Numerical Notes for the definition.
        (
            "EFTpiD1",
            c_double,
        ),  # the value of the pi field equation coefficient D1. Scale independent part of D. Refer to the Numerical Notes for the definition.
        (
            "EFTpiD2",
            c_double,
        ),  # the value of the pi field equation coefficient D2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        (
            "EFTpiE",
            c_double,
        ),  # the value of the pi field equation coefficient E. Refer to the Numerical Notes for the definition.
        # 9) pi field quantities:
        ("pi", c_double),  # the value of the pi field at a given time and scale.
        ("pidot", c_double),  # the value of the (conformal) time derivative of the pi field at a given time and scale.
        (
            "pidotdot",
            c_double,
        ),  # the value of the (conformal) second time derivative of the pi field at a given time and scale.
        # 9b) delta_phi field quantities:
        ("delta_phi", c_double),  # the value of the scalar field perturbation delta_phi
        (
            "delta_phidot",
            c_double,
        ),  # the value of the derivative wrt conformal time of the scalar field perturbation delta_phi
        (
            "delta_phidotdot",
            c_double,
        ),  # the value of the second derivative wrt conformal time of the scalar field perturbation delta_phi
        # 10) scalar perturbations quantities:
        ("etak", c_double),  # Syncronous gauge \f$ eta*k \f$ perturbation.
        (
            "etakdot",
            c_double,
        ),  # Syncronous gauge \f$ \dot{eta}*k \f$ perturbation. This is the Einstein equation that is actually integrated.
        ("z", c_double),  # Syncronous gauge Z perturbation.
        ("dz", c_double),  # Syncronous gauge dot Z perturbation.
        ("sigma", c_double),  # Syncronous gauge sigma perturbation.
        ("sigmadot", c_double),  # Syncronous gauge dot sigma perturbation.
        ("sigmadotdot", c_double),  # Syncronous gauge second time derivative of sigma perturbation.
        ("clxc", c_double),  # Syncronous gauge cdm density perturbation.
        ("clxb", c_double),  # Syncronous gauge baryon density perturbation.
        ("clxg", c_double),  # Syncronous gauge radiation density perturbation.
        ("clxr", c_double),  # Syncronous gauge massless neutrinos density perturbation.
        ("vb", c_double),  # Syncronous gauge baryon velocity.
        ("dgpnu", c_double),  # Syncronous gauge massive neutrinos pressure perturbation.
        ("dgrho", c_double),  # Syncronous gauge total matter density perturbation.
        ("dgq", c_double),  # Syncronous gauge total matter velocity perturbation.
        # 11) tensor perturbations quantities:
        (
            "EFTAT",
            c_double,
        ),  # the value of the tensor equation coefficient A. Refer to the Numerical Notes for the definition.
        (
            "EFTBT",
            c_double,
        ),  # the value of the tensor equation coefficient B. Refer to the Numerical Notes for the definition.
        (
            "EFTDT",
            c_double,
        ),  # the value of the tensor equation coefficient D. Refer to the Numerical Notes for the definition.
        # 12) Kinetic and Gradient quantities for the stability check:
        ("EFT_kinetic", c_double),  # the value of the kinetic term. Refer to the Numerical Notes for the definition.
        ("EFT_gradient", c_double),  # the value of the gradient term. Refer to the Numerical Notes for the definition.
        (
            "EFT_mu1",
            c_double,
        ),  # the value of the first mass eigenvalue. Refer to the Numerical Notes for the definition.
        (
            "EFT_mu2",
            c_double,
        ),  # the value of the second mass eigenvalue. Refer to the Numerical Notes for the definition.
        # 13) other quantities usefull for debug purposes:
        ("Psi", c_double),  # Perturbation in the 00 component of the metric in Conformal Newtonian gauge.
        ("Phi", c_double),  # Perturbation in the space component of the metric in Conformal Newtonian gauge.
        ("PsiDot", c_double),  # Time derivatives of the Newtonian gauge metric perturbations.
        ("PhiDot", c_double),  # Time derivatives of the Newtonian gauge metric perturbations.
        ("EFTISW", c_double),  # Source for ISW effect.
        ("EFTLensing", c_double),  # Source for lensing effect.
        ("T_monopole", c_double),  # The effective CMB temperature monopole.
        ("T_dipole", c_double),  # The effective CMB temperature dipole.
        ("mu", c_double),  # Effective perturbation gravitational constant. \f$ k^2\Psi/\Delta_{m} \f$
        ("gamma", c_double),  # Ratio between the two gravitational potentials. \f$ \gamma = \Phi/\Psi \f$
        (
            "sigma_eft",
            c_double,
        ),  # Ratio between Weyl potential and matter density. \f$ \Sigma = \frac{\Phi+\Psi}{\Delta_{m}} \f$
        (
            "QS_mu",
            c_double,
        ),  # Effective perturbation gravitational constant in the quasi-static limit. \f$ k^2\Psi/\Delta_{m} \f$
        (
            "QS_sigma",
            c_double,
        ),  # Ratio between Weyl potential and matter density in the quasi-static limit. \f$ \Sigma = \frac{\Phi+\Psi}{\Delta_{m}} \f$
        (
            "QS_mu_gen",
            c_double,
        ),  # Effective perturbation gravitational constant in the quasi-static limit. \f$ k^2\Psi/\Delta_{m} \f$
        (
            "QS_sigma_gen",
            c_double,
        ),  # Ratio between Weyl potential and matter density in the quasi-static limit. \f$ \Sigma = \frac{\Phi+\Psi}{\Delta_{m}} \f$
        ("dgrho_Newt", c_double),  # Newtonian gauge total matter density perturbation.
        ("delta_DE_eff", c_double),  # Syncronous gauge Dark energy fluid perturbation.
        # 14) LHS of positivity bounds:
        (
            "posbound1",
            c_double,
        ),  # < LHS of the first positivity bound. Refer to master thesis of Dani de Boe for definition.
        (
            "posbound2",
            c_double,
        ),  # < LHS of the second positivity bound. Refer to master thesis of Dani de Boe for definition.
        (
            "ghostpos",
            c_double,
        ),  # < Ghost condition for the positivity bounds. Refer to master thesis of Dani de Boe for def.
        (
            "tachpos",
            c_double,
        ),  # < Tachyonic condition for the positivity bounds. Refer to master thesis of Dani de Boe for def.
    )

    def dict(self):
        """
        Utility function to return a dictionary of timestep values.
        """
        return {field: getattr(self, field) for field, _ in self._fields_}

    def toarray(self):
        """
        Utility function to return a numpy array with the fields.
        """
        return np.ctypeslib.as_array(self)

    def variables(self):
        """
        Utility to return the names of the fields.
        """
        return [field for field, _ in self._fields_]
