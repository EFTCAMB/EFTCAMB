"""Fast checks of the EFTCAMB Python interface (python -m unittest camb.tests.eftcamb_test).

The EFTCAMB test models are read from fortran/eftcamb_test/parameters (relative to the package) with reduced
accuracy settings (l_max_scalar = 600, no tensors, no source windows, k_max = 2/Mpc), so that each model runs in about
a second. The tests are skipped when the parameter files are not available (e.g. installed wheel).
"""

import copy
import os
import pickle
import re
import shutil
import sys
import tempfile
import unittest
from ctypes import POINTER, byref, c_double, c_int

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb
from camb import eftcamb
from camb.baseconfig import CAMBError, CAMBValueError

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PARAMS = os.path.join(ROOT, "fortran", "eftcamb_test", "parameters")
CACHE_SOURCE = os.path.join(ROOT, "fortran", "eftcamb", "03_abstract_EFTCAMB_cache.f90")

FAST = """
l_max_scalar = 600
get_tensor_cls = F
num_redshiftwindows = 0
DoRedshiftLensing = F
do_lensing = T
lens_potential_accuracy = 0
get_transfer = T
transfer_kmax = 2
transfer_num_redshifts = 1
transfer_redshift(1) = 0
EFTCAMB_write_background = F
feedback_level = 0
"""

# sigma8 and lensed TT C_l (muK^2, raw) at l = 10, 100, 400 with the FAST settings, computed on eftcamb-audit
# (after the audit fixes); relative tolerance 1e-5
REFERENCE = {
    "2_PEFT_Omega_const_1": (0.22928846758790644, [1.52493164e02, 1.65988613e00, 6.38823379e-02]),
    "4_Des_F_R_LCDM_1": (0.9795731927760171, [42.10276032, 1.61639925, 0.06437505]),
    "3_RPH_KT_power_1": (0.6443599338148509, [46.04550157, 1.61509824, 0.06420765]),
    "5_quint_exponential": (0.6283036331682146, [47.14345523, 1.6335124, 0.06406843]),
}


@unittest.skipUnless(os.path.isdir(PARAMS), "EFTCAMB test parameter files not available")
class EFTCAMBTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="eftcamb_test_")
        camb.set_feedback_level(0)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def fast_ini(self, model, extra=""):
        name = os.path.join(self.tmp, f"{model}_{abs(hash(extra))}.ini")
        with open(name, "w") as f:
            f.write(f"DEFAULT({os.path.join(PARAMS, model + '.ini')})\n{FAST}output_root = {self.tmp}/out\n{extra}\n")
        return name

    def read(self, model, extra=""):
        return camb.read_ini(self.fast_ini(model, extra))

    # ctypes layouts

    def test_timestep_cache_layout(self):
        # python must map the Fortran components in the same order
        fields = eftcamb.EFTCAMB_timestep_cache.variables(eftcamb.EFTCAMB_timestep_cache())
        if os.path.exists(CACHE_SOURCE):
            with open(CACHE_SOURCE) as f:
                text = f.read()
            for tname, pyclass in (
                ("TEFTCAMB_timestep_cache", eftcamb.EFTCAMB_timestep_cache),
                ("TEFTCAMB_parameter_cache", eftcamb.EFTCAMB_parameter_cache),
            ):
                body = re.search(
                    rf"\n\s*type\b[^\n]*::\s*{tname}\b(.*?)\n\s*contains", text, re.DOTALL | re.IGNORECASE
                ).group(1)
                names = []
                for line in body.splitlines():
                    code = line.split("!")[0]
                    if "::" in code:
                        names += [re.split(r"[\s=(]", n.strip())[0] for n in code.split("::")[1].split(",")]
                self.assertEqual(
                    [n.lower() for n in names if n], [f.lstrip("_").lower() for f, _ in pyclass._fields_], tname
                )
        # the Fortran code fills an array of timestep caches: with a mismatched element size the fields of later
        # elements would be shifted and the guard element after the array overwritten
        pars = self.read("2_PEFT_Omega_const_1")
        res = camb.get_results(pars)
        a = np.array([0.1, 0.3, 0.6, 1.0])
        eta = np.asfortranarray(res.conformal_time(1 / a - 1))
        n = len(a)
        buf = (eftcamb.EFTCAMB_timestep_cache * (n + 1))()
        guard = np.ctypeslib.as_array(buf)[n : n + 1].view(np.float64)
        guard[:] = 12345.0
        times = np.asfortranarray(a)
        eftcamb.EFTCAMB_eftfunctions(
            byref(res),
            byref(c_int(n)),
            times.ctypes.data_as(POINTER(c_double)),
            eta.ctypes.data_as(POINTER(c_double)),
            byref(buf),
        )
        out = np.ctypeslib.as_array(buf)
        np.testing.assert_allclose(out["a"][:n], a, rtol=1e-14)
        np.testing.assert_allclose(out["tau"][:n], eta, rtol=1e-12)
        np.testing.assert_allclose(out["adotoa"][:n], a * res.hubble_parameter(1 / a - 1) / 299792.458, rtol=1e-6)
        self.assertTrue(np.all(guard == 12345.0))
        self.assertEqual(fields[0], "a")

    def test_parameter_cache_layout(self):
        res = camb.get_results(self.read("2_PEFT_Omega_const_1"))
        p = res.Params
        cache = p.EFTCAMB_parameter_cache
        h = p.H0 / 100
        self.assertAlmostEqual(cache.h0, p.H0, places=10)
        self.assertAlmostEqual(cache.h0_Mpc, p.H0 / 299792.458, places=14)
        self.assertAlmostEqual(cache.omegac, p.omch2 / h**2, places=12)
        self.assertAlmostEqual(cache.omegab, p.ombh2 / h**2, places=12)
        self.assertEqual(cache.Num_Nu_Massive, p.num_nu_massive)
        self.assertEqual(cache.Nu_mass_eigenstates, p.nu_mass_eigenstates)
        self.assertGreaterEqual(len(cache.nu_masses), p.nu_mass_eigenstates)
        self.assertGreater(cache.nu_masses[0], 0)

    def test_parameter_cache_initialised(self):
        # the cache of input parameters is zero (arrays unallocated) until a calculation fills it
        fresh = camb.CAMBparams()
        new = eftcamb.EFTCAMB_parameter_cache()
        read = self.read("2_PEFT_Omega_const_1")
        caches = {"CAMBparams": fresh.EFTCAMB_parameter_cache, "new": new, "read_ini": read.EFTCAMB_parameter_cache}
        for name, cache in caches.items():
            for field, *_ in cache._fields_:
                if field.startswith("_"):
                    self.assertEqual(len(getattr(cache, field[1:])), 0, msg=f"{name}.{field}")
                else:
                    self.assertEqual(getattr(cache, field), 0, msg=f"{name}.{field}")

    # ini reading and the set_params path

    def test_read_ini_set_params(self):
        for model in ("2_PEFT_Omega_power_1", "4_Des_F_R_LCDM_1", "5_quint_exponential"):
            with self.subTest(model=model):
                cp = self.read(model)
                eft = cp.EFTCAMB.read_parameters()
                self.assertNotEqual(eft["EFTflag"], 0)
                self.assertTrue(all(isinstance(v, (bool, int, float, str)) for v in eft.values()))
                cosmo = {"H0": 67.0, "ombh2": 0.022, "omch2": 0.12, "mnu": 0.06}
                a = camb.set_params(**cosmo, **eft)
                b = camb.set_params(**cosmo)
                b.EFTCAMB.initialize_parameters(b, eft, print_header=False)
                self.assertEqual(a.EFTCAMB.model_name(), cp.EFTCAMB.model_name())
                self.assertEqual(a.EFTCAMB.param_values(), cp.EFTCAMB.param_values())
                self.assertEqual(a.EFTCAMB.read_parameters(), b.EFTCAMB.read_parameters())
                z = np.array([0.0, 0.5, 2.0, 10.0])
                np.testing.assert_array_equal(
                    camb.get_background(a).hubble_parameter(z), camb.get_background(b).hubble_parameter(z)
                )
                # reusing the parameters keeps the EFTCAMB model; an EFT-only update changes it
                c = camb.set_params(cp=a, H0=68.0)
                self.assertEqual(c.EFTCAMB.read_parameters(), eft)
                if model == "2_PEFT_Omega_power_1":
                    c = camb.set_params(cp=a, EFTOmega0=0.02)
                    self.assertEqual(c.EFTCAMB.read_parameters()["EFTOmega0"], 0.02)
                    self.assertEqual(c.EFTCAMB.EFTflag, eft["EFTflag"])
        # EFTflag = 0 accepts (and ignores) EFTCAMB keys
        p = camb.set_params(H0=67.0, EFTflag=0, EFT_additional_priors=True)
        self.assertEqual(p.EFTCAMB.EFTflag, 0)

    def test_pickle_deepcopy(self):
        cp = self.read("4_Des_F_R_LCDM_1")
        ref = camb.get_background(cp).hubble_parameter([0.5, 1.0])
        for clone in (pickle.loads(pickle.dumps(cp)), copy.deepcopy(cp), cp.copy()):
            self.assertEqual(clone.EFTCAMB.model_name(), "Designer f(R)")
            self.assertEqual(clone.EFTCAMB.read_parameters(), cp.EFTCAMB.read_parameters())
            np.testing.assert_array_equal(camb.get_background(clone).hubble_parameter([0.5, 1.0]), ref)
        # results.Params can be reused (get_BAO calls the background on a copy)
        res = camb.get_background(cp)
        np.testing.assert_array_equal(res.get_BAO([0.5, 1.0], cp)[:, 1], ref)

    def test_write_ini_round_trip(self):
        cp = self.read("2_PEFT_Omega_power_1")
        name = os.path.join(self.tmp, "written.ini")
        cp.write_ini(name)
        cp2 = camb.read_ini(name)
        self.assertEqual(cp2.EFTCAMB.read_parameters(), cp.EFTCAMB.read_parameters())
        z = np.array([0.0, 1.0, 5.0])
        np.testing.assert_array_equal(
            camb.get_background(cp2).hubble_parameter(z), camb.get_background(cp).hubble_parameter(z)
        )

    # evolution and EFT functions

    def test_get_evolution_order(self):
        cp = self.read("2_PEFT_Omega_const_1")
        res = camb.get_results(cp)
        a = np.array([1.0, 0.05, 0.5, 0.2])
        _, ev = cp.EFTCAMB.get_scale_evolution(res, [0.01, 0.1], a)
        _, ev_sorted = cp.EFTCAMB.get_scale_evolution(res, [0.01, 0.1], np.sort(a))
        self.assertEqual(ev.shape, (2, 4))
        np.testing.assert_allclose(ev["a"][0], a, rtol=1e-8)
        order = np.argsort(a)
        for field in ("a", "tau", "adotoa", "pi"):
            np.testing.assert_array_equal(ev[field][:, order], ev_sorted[field])

    def test_parametrization_derivatives(self):
        # analytic derivatives returned by get_eft_functions vs finite differences of the lower ones
        checks = {
            "2_PEFT_Omega_power_1": [
                ("EFTOmegaV", "EFTOmegaP"),
                ("EFTOmegaP", "EFTOmegaPP"),
                ("EFTOmegaPP", "EFTOmegaPPP"),
            ],
            "2_PEFT_gamma3_exp_1": [("EFTGamma3V", "EFTGamma3P"), ("EFTGamma3P", "EFTGamma3PP")],
        }
        a = np.linspace(0.2, 0.95, 401)
        for model, pairs in checks.items():
            cp = self.read(model)
            _, f = cp.EFTCAMB.get_eft_functions(camb.get_results(cp), a)
            for low, high in pairs:
                with self.subTest(model=model, field=high):
                    fd = np.gradient(f[low], a)[2:-2]
                    np.testing.assert_allclose(
                        fd, f[high][2:-2], rtol=1e-4, atol=1e-6 * np.max(np.abs(f[high])) + 1e-10
                    )
        # Hubble derivatives of a designer and a full-mapping model: dX/dtau = a adotoa dX/da
        for model in ("4_Des_F_R_LCDM_1", "5_quint_exponential"):
            cp = self.read(model)
            _, f = cp.EFTCAMB.get_eft_functions(camb.get_results(cp), a)
            for low, high in (("adotoa", "Hdot"), ("Hdot", "Hdotdot"), ("Hdotdot", "Hdotdotdot")):
                with self.subTest(model=model, field=high):
                    fd = (a * f["adotoa"] * np.gradient(f[low], a))[2:-2]
                    np.testing.assert_allclose(fd, f[high][2:-2], rtol=2e-3, atol=1e-4 * np.max(np.abs(f[high])))

    # errors raise instead of stopping the interpreter

    def test_errors(self):
        cp = camb.set_params(H0=67.0)
        with self.assertRaises(CAMBValueError):
            cp.EFTCAMB.initialize_parameters(cp, {"EFTflag": 1.0}, print_header=False)
        self.assertEqual(cp.EFTCAMB.EFTflag, 0)
        with self.assertRaises(CAMBValueError):
            cp.EFTCAMB.initialize_parameters(
                cp, {"EFTflag": 1, "PureEFTmodel": 1, "EFTCAMB_turn_on_time": -0.1}, print_header=False
            )
        with self.assertRaises(CAMBValueError):
            camb.read_ini(self.fast_ini("2_PEFT_Omega_const_1", "EFTOmega0 = 0.0l"))
        with self.assertRaises(CAMBError):
            camb.get_results(self.read("5_Kmouflage_1", "massive_neutrinos = 1\nomnuh2 = 0.0006"))
        with self.assertRaises(CAMBError):
            camb.get_results(self.read("2_PEFT_Omega_const_1", "omk = -0.05"))
        with self.assertRaises(CAMBError):
            camb.get_results(self.read("2_PEFT_Omega_const_1", "EFTCAMB_evolve_metric_h = T"))
        cp = self.read("2_PEFT_wDE_wCDM_1")
        with self.assertRaises(CAMBValueError):
            cp.set_dark_energy(w=-0.7, dark_energy_model="ppf")
        # GR results have no EFTCAMB model to evaluate
        res = camb.get_background(camb.set_params(H0=67.0))
        with self.assertRaises(CAMBValueError):
            res.Params.EFTCAMB.get_eft_functions(res, [0.5, 1.0])
        # the interpreter is still usable
        self.assertGreater(camb.get_results(self.read("2_PEFT_Omega_const_1")).get_sigma8_0(), 0)

    # GR limit and regression

    def test_gr_flag_zero(self):
        cosmo = {"H0": 67.5, "ombh2": 0.022, "omch2": 0.122, "As": 2e-9, "ns": 0.965, "lmax": 600}
        plain = camb.set_params(**cosmo)
        gr = camb.set_params(**cosmo, EFTflag=0)
        self.assertEqual(gr.EFTCAMB.EFTflag, 0)
        cl_plain = camb.get_results(plain).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:]
        cl_gr = camb.get_results(gr).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:]
        np.testing.assert_allclose(cl_gr[:, :3], cl_plain[:, :3], rtol=1e-10)
        # the shipped GR test model runs through EFTCAMB with EFTflag = 0
        res = camb.get_results(self.read("1_EFT_GR"))
        self.assertEqual(res.Params.EFTCAMB.EFTflag, 0)

    def test_regression(self):
        for model, (sigma8, tt) in REFERENCE.items():
            with self.subTest(model=model):
                res = camb.get_results(self.read(model))
                cl = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][[10, 100, 400], 0]
                self.assertAlmostEqual(res.get_sigma8_0() / sigma8, 1, delta=1e-5)
                np.testing.assert_allclose(cl, tt, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
