"""Fast checks of the EFTCAMB Python interface (python -m unittest camb.tests.eftcamb_test).

The EFTCAMB test models are read from fortran/eftcamb_test/parameters (relative to the package) with reduced
accuracy settings (l_max_scalar = 600, no tensors, no source windows, k_max = 2/Mpc), so that each model runs in about
a second. The tests are skipped when the parameter files are not available (e.g. installed wheel).
"""

import copy
import glob
import os
import pickle
import re
import shutil
import subprocess
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
# (after the audit fixes; the Pure EFT and RPH entries after the pi' initial condition with E', RPH also after the
# 60000-point tables); relative tolerance 1e-5
REFERENCE = {
    "2_PEFT_Omega_const_1": (0.22928843949170907, [1.510155251e02, 1.659871241e00, 6.388233185e-02]),
    "4_Des_F_R_LCDM_1": (0.9795731927760171, [42.10276032, 1.61639925, 0.06437505]),
    "3_RPH_KT_lin_1": (0.6443587523382762, [46.042515957, 1.615098353, 0.064207618]),
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

    def test_qs_outputs_general_form(self):
        # QS_mu and QS_sigma use the general quasi-static form (the former closed form was 12% off for f(R), and
        # ignored gamma3): they equal QS_mu_gen and QS_sigma_gen
        for model in ("4_Des_F_R_LCDM_1", "2_PEFT_gamma3_lin_1"):
            with self.subTest(model=model):
                cp = self.read(model)
                res = camb.get_results(cp)
                _, ev = cp.EFTCAMB.get_scale_evolution(res, [0.005, 0.5, 2.0], [0.1, 0.5, 1.0])
                np.testing.assert_array_equal(ev["QS_mu"], ev["QS_mu_gen"])
                np.testing.assert_array_equal(ev["QS_sigma"], ev["QS_sigma_gen"])
                self.assertTrue(np.all(ev["QS_mu"] != 1))

    def test_shiftsym_m(self):
        # alpha_B = alpha_B0 (H0/H)^(4/m) needs m > 0: Shift_Symmetric_m defaults to 0, and omitting it gave NaN EFT
        # functions and 'theory unstable' with no word on m, even at the GR point alpha_B0 = 0, which now runs as LCDM
        model = "EFTflag = 2\nAltParEFTmodel = 4\nEFTwDE = 0\nShift_Symmetric_alphaK0 = 1.0\n"
        with self.assertRaisesRegex(CAMBError, "Shift_Symmetric_m must be > 0"):
            camb.get_background(self.read("2_PEFT_wDE_LCDM_1", model + "Shift_Symmetric_alphaB0 = 0.1"))
        gr = camb.get_results(self.read("2_PEFT_wDE_LCDM_1", model + "Shift_Symmetric_alphaB0 = 0.0"))
        lcdm = camb.get_results(self.read("2_PEFT_wDE_LCDM_1"))
        self.assertAlmostEqual(gr.get_sigma8_0() / lcdm.get_sigma8_0(), 1, delta=1e-8)

    def test_shiftsym_derivatives(self):
        # shift-symmetric alpha_B: EFTGamma2PP, EFTGamma2PPP and EFTLambdadotdot match finite differences of the
        # lower derivatives (the cache holds a^2 Lambda and a^2 times its conformal-time derivatives)
        model = "EFTflag = 2\nAltParEFTmodel = 4\nShift_Symmetric_alphaK0 = 1.0\n"
        cp = self.read(
            "2_PEFT_wDE_wCDM_1",
            model + "EFTwDE = 1\nEFTw0 = -0.95\nShift_Symmetric_alphaB0 = 0.4\nShift_Symmetric_m = 1.0",
        )
        res = camb.get_results(cp)
        for a0 in (0.3, 0.7, 0.95):
            h = 1e-4 * a0
            e = cp.EFTCAMB.get_scale_evolution(res, [0.1], [a0 - h, a0, a0 + h])[1][0]
            H = e["adotoa"][1]
            fd_g2pp = (e["EFTGamma2P"][2] - e["EFTGamma2P"][0]) / (2 * h)
            fd_g2ppp = (e["EFTGamma2PP"][2] - e["EFTGamma2PP"][0]) / (2 * h)
            fd_ldd = a0 * H * (e["EFTLambdadot"][2] - e["EFTLambdadot"][0]) / (2 * h) - 2 * H * e["EFTLambdadot"][1]
            self.assertLess(abs(e["EFTGamma2PP"][1] - fd_g2pp), 1e-5 * abs(fd_g2pp))
            self.assertLess(abs(e["EFTGamma2PPP"][1] - fd_g2ppp), 1e-5 * abs(fd_g2ppp))
            self.assertLess(abs(e["EFTLambdadotdot"][1] - fd_ldd), 1e-6 * abs(e["EFTLambda"][1]) * H**2)
        # with the exact derivatives the healthy m = 2 model passes the positivity bounds (it used to be refused)
        cp = self.read(
            "2_PEFT_wDE_wCDM_1",
            model + "EFTwDE = 0\nShift_Symmetric_alphaB0 = 0.1\nShift_Symmetric_m = 2.0\nEFT_positivity_bounds = T",
        )
        self.assertGreater(camb.get_results(cp).get_sigma8_0(), 0)

    def test_shiftsym_param_index(self):
        # the w_DE parameters follow alpha_B0, m, alpha_K0: the getters used an offset of 1 and raised
        # 'illegal parameter index' (or returned the wrong entry) for EFTwDE >= 1
        model = "EFTflag = 2\nAltParEFTmodel = 4\nShift_Symmetric_alphaK0 = 1.0\nShift_Symmetric_m = 2.0\n"
        cp = self.read(
            "2_PEFT_wDE_wCDM_1", model + "EFTwDE = 2\nEFTw0 = -0.9\nEFTwa = 0.1\nShift_Symmetric_alphaB0 = 0.1"
        )
        self.assertEqual(cp.EFTCAMB.param_names(), ["alpha_B0", "m", "alpha_K0", "EFTw0", "EFTwa"])
        self.assertEqual(cp.EFTCAMB.param_labels(), ["alpha_B0", "m", "alpha_K0", "w_0", "w_a"])
        self.assertEqual(cp.EFTCAMB.param_values(), [0.1, 2.0, 1.0, -0.9, 0.1])

    def test_braiding_crossing_refused(self):
        # an alpha_B = 2 crossing (Bellini-Sawicki convention) after the turn-on leaves h' undetermined by the 00
        # constraint: refused with its own message (it used to fail with RungeKuttaDP45 -3 or stall)
        model = (
            "EFTflag = 2\nAltParEFTmodel = 4\nShift_Symmetric_alphaK0 = 1.0\nShift_Symmetric_m = 1.0\n"
            "EFTwDE = 1\nEFTw0 = -0.95\n"
        )
        for alphaB0 in ("2.0", "2.1"):
            with self.subTest(alphaB0=alphaB0), self.assertRaisesRegex(CAMBError, "alpha_B = 2 crossing"):
                camb.get_results(self.read("2_PEFT_wDE_wCDM_1", model + f"Shift_Symmetric_alphaB0 = {alphaB0}"))
        cp = self.read("2_PEFT_wDE_wCDM_1", model + "Shift_Symmetric_alphaB0 = 1.98")
        self.assertGreater(camb.get_results(cp).get_sigma8_0(), 0)

    def test_wde_touch_refused(self):
        # quintessence with w_DE touching -1 at one time after the turn-on: the pi kinetic term A = c vanishes there.
        # The touch falls between the sampled scale factors: it is refused cleanly (it used to crash the evolution),
        # while near-touch models and the exact LCDM limit run
        touch = "EFTwa = 0.3\nEFTwat = 0.8\nEFTw0 = "
        with self.assertRaisesRegex(CAMBError, r"vanishes \(A = 0\) at a = 8\.000E-01"):
            camb.get_results(self.read("4_Des_mc5e_turning_point_1", touch + "-1.0"))
        for model, extra in (
            ("4_Des_mc5e_turning_point_1", touch + "-0.999999999999"),
            ("4_Des_mc5e_wCDM_1", "EFTw0 = -1.0"),
        ):
            with self.subTest(model=model, extra=extra):
                self.assertGreater(camb.get_results(self.read(model, extra)).get_sigma8_0(), 0)

    def test_turn_on_warning(self):
        # at feedback_level >= 1 EFTCAMB warns when the requested turn-on cuts into a background on which the GR
        # perturbation equations used before it are inconsistent (|B2|/H^2 > 1e-4); the turn-on stays the user's choice
        code = (
            "import camb, sys\n"
            "for model, extra in (('5_Kmouflage_2', ''), ('5_Kmouflage_2', 'EFTCAMB_turn_on_time = 1e-8'),\n"
            "                     ('2_PEFT_Omega_const_1', '')):\n"
            f"    ini = {self.tmp!r} + '/warn_' + model + str(len(extra)) + '.ini'\n"
            "    open(ini, 'w').write('DEFAULT(' + " + repr(PARAMS) + " + '/' + model + '.ini)\\n' + extra + '\\n')\n"
            "    cp = camb.read_ini(ini)\n"
            "    print('MODEL', model, extra, flush=True)\n"
            "    camb.get_background(cp)\n"
            "    sys.stdout.flush()\n"
        )
        r = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
            env=dict(os.environ, PYTHONPATH=ROOT, OMP_NUM_THREADS="2"),
        )
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        blocks = r.stdout.split("MODEL ")[1:]
        self.assertEqual(len(blocks), 3)
        warning = "the background is not a GR one at the requested EFTCAMB_turn_on_time"
        self.assertIn(warning, blocks[0])
        self.assertNotIn(warning, blocks[1])
        self.assertNotIn(warning, blocks[2])

    def test_glpv_dust_dressed_gradient(self):
        # GLPV gradient condition with pressureless dust dressed (D3): G + (2 G5 - F1)^2 rho_dust/F2^2 >= 0. Pure EFT
        # gamma1 = 0.05, gamma5 = 0.1 checked from a = 0.65: the bare c_s^2 is negative until a ~ 0.7 (-0.23 at 0.66),
        # the dressed one positive (0.34), and the evolved pi propagates with it (c^2 = 0.381 vs 0.384): it was refused.
        # From a = 0.5 the dressed c_s^2 is negative too (-0.32 at 0.56): still refused.
        extra = (
            "PureEFTmodelGamma1 = 1\nEFTGamma10 = 0.05\nPureEFTmodelGamma5 = 1\nEFTGamma50 = 0.1\n"
            "EFT_ghost_math_stability = F\nEFT_ghost_stability = T\nEFT_gradient_stability = T\nEFTCAMB_skip_RGR = T\n"
            "EFTCAMB_turn_on_time = 0.65\nEFTCAMB_stability_time = 0.65"
        )
        self.assertGreater(camb.get_results(self.read("2_PEFT_gamma5_exp_2", extra)).get_sigma8_0(), 0)
        with self.assertRaisesRegex(CAMBError, "theory unstable"):
            camb.get_results(self.read("2_PEFT_gamma5_exp_2", extra.replace("0.65", "0.5")))

    def test_mass_math_deprecated(self):
        # EFT_mass_math_stability is kept as it is (it rejects healthy models: K-mimic here) and warns at feedback >= 1
        code = (
            "import camb, sys\n"
            "for extra in ('EFT_mass_math_stability = T', 'EFT_mass_math_stability = F'):\n"
            f"    ini = {self.tmp!r} + '/mathmass_' + extra[-1] + '.ini'\n"
            "    open(ini, 'w').write('DEFAULT(' + " + repr(PARAMS) + " + '/5_Kmimic_1.ini)\\n' + extra + '\\n')\n"
            "    print('CASE', extra, flush=True)\n"
            "    try:\n"
            "        camb.get_background(camb.read_ini(ini))\n"
            "        print('ACCEPTED', flush=True)\n"
            "    except camb.CAMBError:\n"
            "        print('REFUSED', flush=True)\n"
            "    sys.stdout.flush()\n"
        )
        r = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
            env=dict(os.environ, PYTHONPATH=ROOT, OMP_NUM_THREADS="2"),
        )
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        on, off = r.stdout.split("CASE ")[1:]
        self.assertIn("EFT_mass_math_stability is deprecated", on)
        self.assertIn("REFUSED", on)
        self.assertNotIn("deprecated", off)
        self.assertIn("ACCEPTED", off)

    def test_3p1_gamma2_cosmic_time(self):
        # 3+1 designer: H0 gamma2 = -d rho1/dt with cosmic time, so the model has no braiding (alpha_B = 0:
        # a (H Omega' + H0 gamma2) = 0) and the near-GR model is close to GR (it had a spurious braiding before)
        model = (
            "EFTflag = 3\nDesignerEFTmodel = 3\nEFTwDE = 0\nDesigner3p1LuminalGW = T\nDesigner3p1modelsigma1 = 1\n"
            "Designer3p1sigma10 = 1.0\nDesigner3p1modelsigma1_ODE = 2\nDesigner3p1sigma1_ODE0 = -1.e-4\n"
            "Designer3p1modelrho1N = 0\nDesigner3p1modelrho2 = 1\nDesigner3p1rho20 = 1.e-4\n"
            "EFT_ghost_stability = F\nEFT_gradient_stability = F"
        )
        cp = self.read("2_PEFT_wDE_wCDM_1", model)
        res = camb.get_results(cp)
        h0 = res.Params.H0 / 299792.458
        e = cp.EFTCAMB.get_scale_evolution(res, [0.05], [0.1, 0.5, 0.9])[1][0]
        braiding = e["a"] * (e["adotoa"] * e["EFTOmegaP"] + h0 * e["EFTGamma2V"])
        parts = np.abs(e["a"] * e["adotoa"] * e["EFTOmegaP"]) + np.abs(e["a"] * h0 * e["EFTGamma2V"])
        self.assertTrue(np.all(parts > 0))
        np.testing.assert_array_less(np.abs(braiding), 1e-10 * parts)
        tt = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:601, 0]
        tt_gr = camb.get_results(self.read("1_EFT_GR")).get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"]
        np.testing.assert_allclose(tt, tt_gr[2:601, 0], rtol=1e-3)

    def test_beyond_horndeski_massive_nu(self):
        # the Beyond Horndeski background includes massive neutrinos (base_params has omnuh2 = 0.01): H(z = 0) = H0 and
        # the model H equals the CAMB H (it was 1.0102 H0, and 5e-2 off at early times)
        model = (
            "EFTflag = 4\nFullMappingEFTmodel = 5\nBeyond_Horndeski_x10 = -0.1\nBeyond_Horndeski_x30 = 1.0\n"
            "Beyond_Horndeski_x40 = 1.0\nEFTCAMB_stability_threshold = 0.0"
        )
        cp = self.read("2_PEFT_wDE_wCDM_1", model)
        self.assertGreater(cp.omnuh2, 0.005)
        res = camb.get_results(cp)
        self.assertAlmostEqual(res.hubble_parameter(0.0) / cp.H0, 1, delta=1e-5)
        a = np.array([1e-4, 1e-2, 0.5, 1.0])
        h_eft = cp.EFTCAMB.get_eft_functions(res, a)[1]["adotoa"]
        np.testing.assert_allclose(h_eft, a * res.h_of_z(1 / a - 1), rtol=1e-5)

    def test_beyond_horndeski_lambda(self):
        # Beyond Horndeski EFTLambda has the 3 x4/5 term (it had 6 x4/5): CAMB's H, from the EFT functions, equals
        # the model H (they differed by -x4/10, 2.6e-4 for the PBFT19 BH1 benchmark)
        for x40 in ("1.0", "100.0"):
            with self.subTest(x40=x40):
                x10 = "-0.1" if x40 == "1.0" else "-1.0"
                model = (
                    f"EFTflag = 4\nFullMappingEFTmodel = 5\nBeyond_Horndeski_x10 = {x10}\nBeyond_Horndeski_x30 = 1.0\n"
                    f"Beyond_Horndeski_x40 = {x40}\nEFTCAMB_stability_threshold = 0.0"
                )
                cp = self.read("2_PEFT_wDE_wCDM_1", model)
                res = camb.get_results(cp)
                a = np.array([1e-4, 1e-2, 0.5, 1.0])
                h_eft = cp.EFTCAMB.get_eft_functions(res, a)[1]["adotoa"]
                np.testing.assert_allclose(h_eft, a * res.h_of_z(1 / a - 1), rtol=1e-12)

    def test_beyond_horndeski_x4_zero_derivatives(self):
        # with x4 = 0 the Beyond Horndeski branch had four typos in the higher derivatives (OmegaRad for OmegaRad', 3 x3 x3'
        # for 3 (x3'^2 + x3 x3''), 9 x2'' for 9 x1'', Hdot^3 x3 for 3 Hdot^3 x3): EFTLambdadotdot was 13-63x off,
        # Hdotdotdot 1-4x, EFTGamma2PPP up to 9%. Against central differences of the next lower derivative (the
        # Lambda fields are a^2-weighted: a^2 Lambda'' = d(a^2 Lambda')/dtau - 2 Hc a^2 Lambda'):
        model = "EFTflag = 4\nFullMappingEFTmodel = 5\nBeyond_Horndeski_x10 = -0.1\nBeyond_Horndeski_x30 = 1.0\n"
        cp = self.read("2_PEFT_wDE_wCDM_1", model + "Beyond_Horndeski_x40 = 0.0\nEFTCAMB_skip_stability = T")
        bg = camb.get_background(cp)
        names = ("adotoa", "Hdotdot", "Hdotdotdot", "EFTLambdadot", "EFTLambdadotdot", "EFTGamma2PP", "EFTGamma2PPP")
        for a0 in (0.3, 0.7, 0.95):
            h = 1e-4 * a0
            f = cp.EFTCAMB.get_eft_functions(bg, np.array([a0 - h, a0, a0 + h]))[1]
            g = {k: np.ravel(f[k]) for k in names}
            dg = {k: (g[k][2] - g[k][0]) / (2 * h) for k in names}
            Hc = g["adotoa"][1]
            with self.subTest(a=a0):
                self.assertAlmostEqual(g["Hdotdotdot"][1] / (a0 * Hc * dg["Hdotdot"]), 1, delta=5e-3)
                lam_fd = a0 * Hc * dg["EFTLambdadot"] - 2 * Hc * g["EFTLambdadot"][1]
                self.assertAlmostEqual(g["EFTLambdadotdot"][1] / lam_fd, 1, delta=5e-3)
                self.assertAlmostEqual(g["EFTGamma2PPP"][1] / dg["EFTGamma2PP"], 1, delta=5e-3)

    def test_beyond_horndeski_background_solver(self):
        # no hidden Omega_v >= 0.3 cut (omch2 = 0.35 was refused as 'background solver failed'); before the background
        # table (a_s = 1/(1.5e5 + 1)) x_i grow with the rates of the solution at a_s: for x4 = 0 the fixed x4-attractor
        # exponents made d ln c/d ln a jump from +1.0 to -1.5 at a_s; a failed integration is named
        model = "EFTflag = 4\nFullMappingEFTmodel = 5\nBeyond_Horndeski_x10 = -0.1\nBeyond_Horndeski_x30 = 1.0\n"
        model += "EFTCAMB_stability_threshold = 0.0\nBeyond_Horndeski_x40 = "
        cp = self.read("2_PEFT_wDE_wCDM_1", model + "1.0\nomch2 = 0.35")
        self.assertAlmostEqual(camb.get_results(cp).hubble_parameter(0.0) / cp.H0, 1, delta=1e-5)
        cp = self.read("2_PEFT_wDE_wCDM_1", model + "0.0")
        bg = camb.get_background(cp)
        a = np.array([0.95, 1.05]) / (1.5e5 + 1)
        h = 0.01
        c = [np.array(cp.EFTCAMB.get_eft_functions(bg, a * np.exp(s * h))[1]["EFTc"]) for s in (-1, 1)]
        slope = np.log(c[1] / c[0]) / (2 * h)
        self.assertAlmostEqual(slope[0], slope[1], delta=1e-2)
        with self.assertRaisesRegex(CAMBError, "Beyond Horndeski background integration failed"):
            camb.get_results(self.read("2_PEFT_wDE_wCDM_1", model + "-100.0"))

    def test_ade_flatness(self):
        # ADE closes flatness with the continuity integral for the input p (the p = 1 closed form gave H(z = 0) 2.5% off
        # for p = 1/2, w_f = 0)
        for extra in (
            "wf = 0.0\nLog_ac = -3.5\nf_ac = 0.1\np = 0.5",
            "wf = 0.1\nLog_ac = -3.5\nf_ac = 0.1\np = 0.5",
            "",
        ):
            with self.subTest(extra=extra):
                cp = self.read("5_ADE_1", extra)
                self.assertAlmostEqual(camb.get_background(cp).hubble_parameter(0.0) / cp.H0, 1, delta=1e-9)

    def test_ade_sound_speed(self):
        # ADE evaluates 1 + w directly: 1 + (w - 1) rounded to exactly 0 deep before a_c, so c = gamma1 = 0 while c' and
        # Lambda' were not, and this c_s^2 != 1 model was refused as having no pi kinetic term
        a = np.logspace(-8, -6, 5)
        cp = self.read("5_ADE_1", "cs2 = 0.86\nwf = 0.87\nLog_ac = -3.52\nf_ac = 0.086\np = 0.5")
        res = camb.get_results(cp)
        ev = cp.EFTCAMB.get_scale_evolution(res, [0.05], a)[1][0]
        c, gamma1 = np.array(ev["EFTc"]), np.array(ev["EFTGamma1V"])
        self.assertTrue(np.all(c > 0))
        h0 = cp.H0 / 299792.458  # H0 in 1/Mpc
        np.testing.assert_allclose(c / (c + 2 * h0**2 * a**2 * gamma1), 0.86, rtol=1e-7)

    def test_ade_priors(self):
        # ADE: the undefined domain (f_ac outside [0, 1), p <= 0) is refused by the background initialization, the
        # priors with named messages (they only printed, and cs2 < 0 ran as 'stable' and failed in the integrator).
        # f_ac = 0 (no ADE) is LCDM and runs: it was refused by the f_ac prior.
        for extra, message in (
            ("cs2 = -0.5", "cs2 < 0"),
            ("p = -1.0", "p <= 0"),
            ("f_ac = 1.0", r"f_ac outside \[0, 1\)"),
            ("Log_ac = 0.5", "Log_ac > 0"),
        ):
            with self.subTest(extra=extra), self.assertRaisesRegex(CAMBError, message):
                camb.get_background(self.read("5_ADE_1", extra))
        lcdm = camb.get_results(self.read("2_PEFT_wDE_LCDM_1")).get_sigma8_0()
        ade = camb.get_results(self.read("5_ADE_1", "f_ac = 0.0")).get_sigma8_0()
        self.assertAlmostEqual(ade / lcdm, 1, delta=1e-8)

    def test_ade_turn_on(self):
        # the return to GR compares d gamma_i/d ln a with the threshold: the raw d gamma1/da = 1e-8 switched this ADE
        # (arXiv:1905.12618 Table I) on at a = 7.3e-9, where gamma1 = 1e-17
        cp = self.read("5_ADE_1", "cs2 = 0.86\nwf = 0.87\nLog_ac = -3.52\nf_ac = 0.086\np = 0.5")
        bg = camb.get_background(cp)
        self.assertGreater(bg.Params.EFTCAMB.EFTCAMB_pert_turn_on, 1e-7)

    def test_ic_pidot_attractor(self):
        # the type-1 initial pi' is the time derivative of the attractor -H0 E/(C + k^2 D1 + k^4 D2) (arXiv:1405.3590),
        # with E' for every model: E' = 0 whenever a gamma_i != 0 gave pi' = -4x (RPH) and ~0 (Horava) that derivative
        eps = np.array([1e-6, 2e-6, 3e-6])
        for model in ("3_RPH_KT_OmegaDE_1", "5_Horava_eta_3"):
            cp = self.read(model)
            res = camb.get_results(cp)
            h0 = cp.H0 / 299792.458  # H0 in 1/Mpc
            tau_on = res.conformal_time(1 / res.Params.EFTCAMB.EFTCAMB_pert_turn_on - 1)
            for k in (0.01, 0.1, 1.0):
                with self.subTest(model=model, k=k):
                    ev = cp.EFTCAMB.get_evolution(res, [k], tau_on * (1 + eps))[1][0]
                    den = ev["EFTpiC"] + k**2 * ev["EFTpiD1"] + k**4 * ev["EFTpiD2"]
                    attractor = -h0 * ev["EFTpiE"] / den
                    dattractor = (attractor[2] - attractor[0]) / (tau_on * (eps[2] - eps[0]))
                    self.assertAlmostEqual(ev["pidot"][1] / dattractor, 1, delta=0.1)

    def test_quintessence_h0_shooting(self):
        # the H0 shooting returns a root that reproduces H0 (with drag, p = 1 returned H0 3.4% off; without drag it was
        # refused), the same model with and without drag, and refuses a steep exponential whose only root starts
        # scalar dominated (Omega_phi(a_ini) = 1: no radiation or matter era)
        sigma8 = []
        for drag in ("T", "F"):
            cp = self.read("5_quint_quadratic", f"p = 1.0\ndrag_initial_conditions = {drag}")
            res = camb.get_results(cp)
            self.assertAlmostEqual(res.hubble_parameter(0) / cp.H0, 1, delta=1e-6)
            sigma8.append(res.get_sigma8_0())
        self.assertAlmostEqual(sigma8[0] / sigma8[1], 1, delta=1e-7)
        with self.assertRaisesRegex(CAMBError, "background solver failed"):
            camb.get_results(self.read("5_quint_exponential", "EFT_lambda = 2.2"))

    def test_evolution_higher_derivatives(self):
        # get_scale_evolution computes the higher derivatives before the mass eigenvalues, as the stability check does:
        # for full models they were 0 there (EFTcdotdot, EFTLambdadotdot) and Hdotdotdot was 3x off at a = 0.5
        a = np.array([0.02, 0.1, 0.5])
        cp = self.read("5_quint_exponential")
        res = camb.get_results(cp)
        ev = cp.EFTCAMB.get_scale_evolution(res, [0.05], a)[1][0]
        eft = cp.EFTCAMB.get_eft_functions(res, a)[1]
        for name in ("EFTcdotdot", "EFTLambdadotdot", "Hdotdotdot"):
            with self.subTest(name=name):
                np.testing.assert_allclose(ev[name], eft[name], rtol=1e-6)

    def test_quintessence_lcdm_point(self):
        # a constant potential with V0 = Omega_v/(1 - Om_c - Om_b - Om_nu) (the potential unit) is LCDM: the a^2-weighted
        # Lambda table switched it on at a = 0.17, where it was refused (no pi kinetic term); it is never switched on now
        lcdm = camb.get_results(self.read("2_PEFT_wDE_LCDM_1"))
        cp = self.read("5_quint_quadratic")
        v0 = lcdm.get_Omega("de") / (1 - (cp.omch2 + cp.ombh2 + cp.omnuh2) / (cp.H0 / 100) ** 2)
        extra = f"potential_model = 1\nV0 = {float(v0)!r}\ndrag_initial_conditions = F\nphidot_ini = 0.0"
        res = camb.get_results(self.read("5_quint_quadratic", extra))
        self.assertGreater(res.Params.EFTCAMB.EFTCAMB_pert_turn_on, 1)
        self.assertAlmostEqual(res.get_sigma8_0() / lcdm.get_sigma8_0(), 1, delta=1e-10)

    def test_quintessence_early_background(self):
        # dphi/dN is 1e-28..1e-9 before a ~ 1e-3: with one absolute tolerance 1e-14 (and a transposed Jacobian) it was
        # solver noise, and EFTcdot/(H EFTc) = 38 at a = 1e-6 while c grew as a^4. The cache holds a^2 c and a^2 c'.
        a = np.array([1e-7, 1e-6, 1e-5, 1e-4])
        h = 0.05
        cp = self.read("5_quint_exponential")
        bg = camb.get_background(cp)
        f = cp.EFTCAMB.get_eft_functions(bg, a)[1]
        slope = np.array(f["EFTcdot"]) / (np.array(f["adotoa"]) * np.array(f["EFTc"]))
        c_plus = np.array(cp.EFTCAMB.get_eft_functions(bg, a * np.exp(h))[1]["EFTc"])
        c_minus = np.array(cp.EFTCAMB.get_eft_functions(bg, a * np.exp(-h))[1]["EFTc"])
        np.testing.assert_allclose(slope, np.log(c_plus / c_minus) / (2 * h) - 2, atol=1e-3)

    def test_scaling_cubic_scaling_point(self):
        # the initial conditions are the radiation scaling point, which exists only if 4/3 + 2 A (beta1 - 4 lambda/3) >= 0:
        # A = 0.1 gave NaN initial conditions and 'Scaling Cubic background is NaN'
        with self.assertRaisesRegex(CAMBError, "radiation scaling point does not exist"):
            camb.get_background(self.read("5_scg", "Scaling_Cubic_A = 0.1"))

    def test_scaling_cubic_shooting(self):
        # the y2_i bracket expands until Omega_DE is overshot, and there is no Omega_v >= 0.3 cut: these were refused
        for extra in ("Scaling_Cubic_beta2 = 3.0", "model_background_a_ini = 1e-5", "omch2 = 0.33", "hubble = 99.0"):
            with self.subTest(extra=extra):
                cp = self.read("5_scg", extra)
                self.assertAlmostEqual(camb.get_results(cp).hubble_parameter(0) / cp.H0, 1, delta=1e-5)

    def test_scaling_cubic_early_background(self):
        # the Scaling Cubic background starts at a = 1e-10: it started at 6.7e-6, after the 1e-8 turn-on of 5_scg, and
        # was frozen there (gamma2' 4.2e-4 off d gamma2/da at a = 1e-6). A turn-on before model_background_a_ini warns.
        a, h = np.array([1e-8, 1e-6]), 1e-3
        cp = self.read("5_scg")
        bg = camb.get_background(cp)
        g2 = [np.array(cp.EFTCAMB.get_eft_functions(bg, a * np.exp(s * h))[1]["EFTGamma2V"]) for s in (-1, 1)]
        g2p = np.array(cp.EFTCAMB.get_eft_functions(bg, a)[1]["EFTGamma2P"])
        np.testing.assert_allclose(g2p, (g2[1] - g2[0]) / (2 * h * a), rtol=1e-5)
        code = (
            "import camb\n"
            "for extra in ('', 'model_background_a_ini = 1e-6'):\n"
            f"    ini = {self.tmp!r} + '/scg_start' + str(len(extra)) + '.ini'\n"
            "    open(ini, 'w').write('DEFAULT(' + "
            + repr(PARAMS)
            + " + '/5_scg.ini)\\nfeedback_level = 1\\n' + extra)\n"
            "    print('MODEL', extra, flush=True)\n"
            "    camb.get_background(camb.read_ini(ini))\n"
        )
        env = dict(os.environ, PYTHONPATH=ROOT, OMP_NUM_THREADS="2")
        r = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, timeout=600, env=env, check=False
        )
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        blocks = r.stdout.split("MODEL ")[1:]
        self.assertEqual(len(blocks), 2)
        warning = "before the model background starts"
        self.assertNotIn(warning, blocks[0])
        self.assertIn(warning, blocks[1])

    def test_kmouflage_parameter_domain(self):
        # K-mimic does not use alphaU (alphaU < 0 was refused); parameters for which the background is undefined (they
        # were refused through a NaN dtauda) are named, also without the priors, and so are the priors
        sigma8 = camb.get_results(self.read("5_Kmimic_1")).get_sigma8_0()
        self.assertEqual(camb.get_results(self.read("5_Kmimic_1", "alphaU = -0.1")).get_sigma8_0(), sigma8)
        for model, extra, message in (
            ("5_Kmouflage_3", "eps2_0 = 0.0\nEFT_additional_priors = F", r"undefined for eps2_0 = 0 \(the LCDM limit"),
            ("5_Kmimic_1", "gammaA = 0.0", r"undefined for gammaA <= 0"),
            ("5_Kmimic_1", "m = 1.0", r"undefined for m <= 1"),
            ("5_Kmimic_1", "m = 1.2\neps2_0 = 0.1\ngammaA = 0.2", r"undefined for A\(a -> 0\) <= 0"),
            ("5_Kmouflage_3", "alphaU = -0.1", r"prior not met: alphaU < 0"),
            ("5_Kmimic_1", "eps2_0 = 0.1\ngammaA = 0.1", r"prior not met: gammaA < eps2_0/\(1 - eps2_0\)"),
        ):
            with self.subTest(extra=extra), self.assertRaisesRegex(CAMBError, message):
                camb.get_results(self.read(model, extra))

    def test_extended_galileon_prior(self):
        # the empirical (S + 0.02) q < 0.5 line (a fit to the optional mass cut of arXiv:1912.07586) refused these
        # viable models; the strong-coupling prior still refuses q S > 1
        model = "EFTflag = 4\nFullMappingEFTmodel = 7\nEFT_ghost_stability = T\nEFT_gradient_stability = T\n"
        model += "EFT_additional_priors = T\nEFTCAMB_stability_threshold = 0.0\n"
        for s, q in (("0.05", "1.0"), ("1.0", "0.35")):
            with self.subTest(S=s, q=q):
                cp = self.read("1_EFT_GR", f"{model}ExtendedGalileon_S = {s}\nExtendedGalileon_q = {q}")
                self.assertGreater(camb.get_results(cp).get_sigma8_0(), 0)
        with self.assertRaises(CAMBError):
            camb.get_results(self.read("1_EFT_GR", f"{model}ExtendedGalileon_S = 2.0\nExtendedGalileon_q = 1.0"))
        # the boundary q S = 1 is allowed: (0.1, 10) and (0.2, 5) were refused by roundoff (q S - 1 = 5.6e-17 with
        # fused multiply-add); q S = 1.1 is still refused
        for s, q in (("0.1", "10.0"), ("0.2", "5.0")):
            with self.subTest(S=s, q=q):
                camb.get_background(self.read("1_EFT_GR", f"{model}ExtendedGalileon_S = {s}\nExtendedGalileon_q = {q}"))
        with self.assertRaisesRegex(CAMBError, "theory unstable"):
            camb.get_background(self.read("1_EFT_GR", f"{model}ExtendedGalileon_S = 0.11\nExtendedGalileon_q = 10.0"))

    def test_extended_galileon_tracker_cache(self):
        # the tracker solution is cached per thread, keyed on every input: another model evaluated just before at
        # the same a (other S; other split of the same omnuh2 between two eigenstates, which a key on the sum of
        # m_nu grhormass would miss; other H0) must not change the result. The Fortran routine is called directly,
        # since get_eft_functions first integrates the conformal time of the model (which refills the cache)
        def eft_at(res):
            buf = (eftcamb.EFTCAMB_timestep_cache * 1)()
            a, eta = np.array([0.5]), np.array([0.0])
            ptr = (a.ctypes.data_as(POINTER(c_double)), eta.ctypes.data_as(POINTER(c_double)))
            eftcamb.EFTCAMB_eftfunctions(byref(res), byref(c_int(1)), *ptr, byref(buf))
            return np.ctypeslib.as_array(buf)

        model = "EFTflag = 4\nFullMappingEFTmodel = 7\nExtendedGalileon_q = 0.35\nmassless_neutrinos = 1.046\n"
        model += "nu_mass_eigenstates = 2\nmassive_neutrinos = 1 1\nExtendedGalileon_S = "
        bg = camb.get_background(self.read("1_EFT_GR", model + "1.0\nnu_mass_fractions = 0.5 0.5"))
        ref = eft_at(bg)
        for extra in (
            "1.5\nnu_mass_fractions = 0.5 0.5",
            "1.0\nnu_mass_fractions = 0.3 0.7",
            "1.0\nnu_mass_fractions = 0.5 0.5\nhubble = 67",
        ):
            with self.subTest(extra=extra):
                other = eft_at(camb.get_background(self.read("1_EFT_GR", model + extra)))
                f = eft_at(bg)
                self.assertNotEqual(other["adotoa"][0], ref["adotoa"][0])
                self.assertEqual(f.tobytes(), ref.tobytes())

    def test_omega_lambda_backgrounds(self):
        # OL gamma / OL alpha integrate the deviation from the LCDM solution of the same matter from a = 1 (no
        # shooting, no dust mode): the GR points (refused before) are the Pure EFT LCDM, and H(z = 0) = H0
        lcdm = camb.get_results(self.read("2_PEFT_wDE_LCDM_1")).get_sigma8_0()
        olg = "EFTflag = 2\nAltParEFTmodel = 2\nOLLambdamodel = 0\nOLGamma1model = 0\nOLGamma2model = 0\n"
        olg += "OLGamma3model = 0\nOLOmegamodel = "
        ola = "EFTflag = 2\nAltParEFTmodel = 3\nOLLambdamodel = 0\nOLmassPmodel = 0\nOLkineticitymodel = 0\n"
        ola += "OLbraidingmodel = 0\nOLtensormodel = 0"
        for extra in (olg + "0", ola):
            with self.subTest(extra=extra):
                res = camb.get_results(self.read("2_PEFT_wDE_LCDM_1", extra))
                self.assertAlmostEqual(res.get_sigma8_0() / lcdm, 1, delta=1e-10)
        cp = self.read("2_PEFT_wDE_LCDM_1", olg + "3\nOLOmega0 = 0.05\nOLOmegaExp = 1.0")
        self.assertAlmostEqual(camb.get_results(cp).hubble_parameter(0) / cp.H0, 1, delta=1e-7)

    def test_omega_lambda_alpha_gamma3_derivatives(self):
        # OL alpha sets the third derivative of gamma3 = -(1 + M) alpha_T; with M = 0.1 a, alpha_T = -0.05 a^3 the third
        # and fourth derivatives of gamma3 = 0.05 a^3 + 0.005 a^4 were both 0
        ola = "EFTflag = 2\nAltParEFTmodel = 3\nOLLambdamodel = 0\nOLmassPmodel = 3\nOLkineticitymodel = 3\n"
        ola += "OLbraidingmodel = 3\nOLtensormodel = 3\nOLmass0 = 0.1\nOLmassExp = 1.0\nOLkineticity0 = 1.0\n"
        ola += "OLkineticityExp = 1.0\nOLbraiding0 = 0.1\nOLbraidingExp = 1.0\nOLtensor0 = -0.05\nOLtensorExp = 3.0\n"
        ola += "EFTCAMB_skip_stability = T"
        cp = self.read("2_PEFT_wDE_LCDM_1", ola)
        a = np.array([0.3, 0.6, 0.9])
        f = cp.EFTCAMB.get_eft_functions(camb.get_background(cp), a)[1]
        np.testing.assert_allclose(f["EFTGamma3V"], 0.05 * a**3 + 0.005 * a**4, rtol=1e-12)
        np.testing.assert_allclose(f["EFTGamma3PPP"], 0.3 + 0.12 * a, rtol=1e-12)
        np.testing.assert_allclose(f["EFTGamma3PPPP"], 0.12, rtol=1e-6)

    def test_numerical_derivative_fallbacks(self):
        # near a = 1 the fallback derivatives are one-sided (the central stencil went past the end of the tables:
        # designer f(R) OmegaPPPP(1) = -0.46 off the trend -0.92); a third derivative of a field whose second
        # derivative the model does not set is the second derivative of the first (EFTGamma2PPP was 0 for OL alpha and
        # Extended Galileon)
        cp = self.read("4_Des_F_R_LCDM_1")
        f = cp.EFTCAMB.get_eft_functions(camb.get_background(cp), np.array([0.97, 0.98, 0.99, 1.0]))[1]
        for name in ("EFTOmegaPPPP", "EFTLambdadotdot"):
            v = np.array(f[name])
            with self.subTest(name=name):
                self.assertAlmostEqual(v[3] / (3 * v[2] - 3 * v[1] + v[0]), 1, delta=1e-2)
        ola = "EFTflag = 2\nAltParEFTmodel = 3\nOLLambdamodel = 0\nOLmassPmodel = 3\nOLkineticitymodel = 3\n"
        ola += "OLbraidingmodel = 3\nOLtensormodel = 3\nOLmass0 = 0.1\nOLmassExp = 1.0\nOLkineticity0 = 1.0\n"
        ola += "OLkineticityExp = 1.0\nOLbraiding0 = 0.1\nOLbraidingExp = 1.0\nOLtensor0 = -0.05\nOLtensorExp = 1.0\n"
        galileon = "EFTflag = 4\nFullMappingEFTmodel = 7\nExtendedGalileon_S = 0.05\nExtendedGalileon_q = 1.0\n"
        a, h = np.array([0.3, 0.6, 0.9]), 1e-3
        for model, extra in (("2_PEFT_wDE_LCDM_1", ola), ("1_EFT_GR", galileon)):
            with self.subTest(model=extra[:30]):
                cp = self.read(model, extra + "EFTCAMB_skip_stability = T")
                bg = camb.get_background(cp)
                pp = [np.array(cp.EFTCAMB.get_eft_functions(bg, a * (1 + s * h))[1]["EFTGamma2PP"]) for s in (-1, 1)]
                ppp = np.array(cp.EFTCAMB.get_eft_functions(bg, a)[1]["EFTGamma2PPP"])
                np.testing.assert_allclose(ppp, (pp[1] - pp[0]) / (2 * h * a), rtol=1e-3)

    def test_numerical_derivative_large_values(self):
        # the numerical derivative (dfridr) of the fallbacks returned 0 for derivatives above ~1e30: OL alpha
        # EFTGamma1PP (-6.7e39 .. -6.9e32 for a = 1e-8 .. 2.5e-7) was 0 there
        ola = "EFTflag = 2\nAltParEFTmodel = 3\nOLLambdamodel = 0\nOLmassPmodel = 3\nOLkineticitymodel = 3\n"
        ola += "OLbraidingmodel = 0\nOLtensormodel = 0\nOLmass0 = 0.05\nOLmassExp = 1.0\nOLkineticity0 = 1.0\n"
        ola += "OLkineticityExp = 1.0\nEFTCAMB_skip_stability = T"
        cp = self.read("2_PEFT_wDE_LCDM_1", ola)
        bg = camb.get_background(cp)
        a, h = np.array([1e-8, 1e-7, 2.5e-7, 1e-6]), 1e-3
        g1p = [np.array(cp.EFTCAMB.get_eft_functions(bg, a * (1 + s * h))[1]["EFTGamma1P"]) for s in (-1, 1)]
        g1pp = np.array(cp.EFTCAMB.get_eft_functions(bg, a)[1]["EFTGamma1PP"])
        np.testing.assert_allclose(g1pp, (g1p[1] - g1p[0]) / (2 * h * a), rtol=1e-4)

    def test_glpv_domain_floor(self):
        # gamma6 and gamma3 + gamma4 enter the dispersion relation times (k/aH)^2 ~ 1e9: only roundoff counts as GLPV
        # (the physical floor 1e-8 let gamma3 + gamma4 = -1e-9, with a k^4 instability above 2.2/Mpc, pass)
        model = "PureEFTHorndeski = F\nPureEFTmodelGamma"
        for extra in ("3 = 1\nEFTGamma30 = -1e-9", "6 = 1\nEFTGamma60 = -1e-9"):
            with self.subTest(extra=extra), self.assertRaisesRegex(CAMBError, "theory unstable"):
                camb.get_results(self.read("2_PEFT_wDE_LCDM_1", model + extra))
        self.assertGreater(
            camb.get_results(self.read("2_PEFT_wDE_LCDM_1", model + "3 = 1\nEFTGamma30 = -1e-13")).get_sigma8_0(), 0
        )

    def test_fR_gr_limit(self):
        # designer f(R) with B0 = 0 is LCDM: the return to GR never switches it on (the interpolation error of an
        # a^2-weighted Lambda table used to switch it on at a = 0.099, then the run failed)
        lcdm = camb.get_results(self.read("2_PEFT_wDE_LCDM_1"))
        for b0 in ("0.0", "1.e-10"):
            with self.subTest(B0=b0):
                res = camb.get_results(self.read("4_Des_F_R_LCDM_1", f"EFTB0 = {b0}"))
                self.assertGreater(res.Params.EFTCAMB.EFTCAMB_pert_turn_on, 1)
                self.assertAlmostEqual(res.get_sigma8_0() / lcdm.get_sigma8_0(), 1, delta=1e-12)

    def test_fR_b0_shooting(self):
        # the B0 shooting must reproduce the requested B0 (relative 1e-4), or the model is refused: w = -0.45, B0 = 1
        # returns B(1) = 1.001. w = -0.55, -0.8 with B0 = 1 (missed while the initial condition counted massive
        # neutrinos as radiation: B(1) = 6.1e5 for -0.55) and w = -0.85, B0 = 1e-4 run
        with self.assertRaisesRegex(CAMBError, "background solver failed"):
            camb.get_results(self.read("4_Des_F_R_wCDM_1", "EFTw0 = -0.45"))
        for extra in ("EFTw0 = -0.55", "EFTw0 = -0.8", "EFTw0 = -0.85\nEFTB0 = 1e-4"):
            with self.subTest(extra=extra):
                self.assertGreater(camb.get_results(self.read("4_Des_F_R_wCDM_1", extra)).get_sigma8_0(), 0)

    def test_fR_massive_nu_ic(self):
        # the f(R) initial condition at a = 1e-8 counts massive neutrinos as radiation, as they are there: Omega(a) then
        # starts as with 3.046 massless neutrinos (today's omegan as radiation gave d ln Omega/d ln a = 958 at a_ini)
        a = np.array([1e-8, 1.2e-8, 1.5e-8, 2e-8])
        slopes = []
        for extra in ("", "omnuh2 = 0\nmassless_neutrinos = 3.046\nmassive_neutrinos = 0"):
            cp = self.read("4_Des_F_R_LCDM_1", extra)
            f = cp.EFTCAMB.get_eft_functions(camb.get_background(cp), a)[1]
            slopes.append(a * np.array(f["EFTOmegaP"]) / np.array(f["EFTOmegaV"]))
        np.testing.assert_allclose(slopes[0], slopes[1], rtol=1e-5)

    def test_fR_tachyon_prior(self):
        # the f(R) prior refuses f_RR < 0 (tachyonic scalaron) from EFTCAMB_stability_time on: B0 < 0 and a constant
        # w < -1 (B < 0 only before a = 3.4e-3 for w = -1.01, B0 = 1: it runs if the check starts later)
        for model, extra in (("4_Des_F_R_LCDM_1", "EFTB0 = -1e-4"), ("4_Des_F_R_wCDM_1", "EFTw0 = -1.01")):
            with self.subTest(extra=extra), self.assertRaisesRegex(CAMBError, r"f_RR < 0 \(tachyonic scalaron\)"):
                camb.get_results(self.read(model, extra))
        res = camb.get_results(self.read("4_Des_F_R_wCDM_1", "EFTw0 = -1.01\nEFTCAMB_stability_time = 1e-2"))
        self.assertGreater(res.get_sigma8_0(), 0)

    def fR_mass(self, extra, a):
        """EFT_mu1 of designer f(R) from get_scale_evolution, and the scalaron mass ((1 + f_R)/f_RR - R)/3 from the same
        cache (f_RR = dOmega/dt / dR/dt, cosmic time)."""
        cp = self.read("4_Des_F_R_LCDM_1", extra)
        ev = cp.EFTCAMB.get_scale_evolution(camb.get_results(cp), [1e-3], a)[1][0]
        ev = {name: np.ravel(ev[name]) for name in ("adotoa", "Hdot", "Hdotdot", "EFTOmegaV", "EFTOmegaP", "EFT_mu1")}
        H = ev["adotoa"] / a
        Hd = ev["Hdot"] / a**2 - H**2
        Hdd = ev["Hdotdot"] / a**3 - 4 * H * Hd - 2 * H**3
        R, Rdot = 6 * (2 * H**2 + Hd), 6 * (4 * H * Hd + Hdd)
        return ev["EFT_mu1"], ((1 + ev["EFTOmegaV"]) * Rdot / (ev["adotoa"] * ev["EFTOmegaP"]) - R) / 3

    def test_alpha_functions(self):
        # get_eft_functions and the evolution fill the alpha-basis outputs from the EFT functions for every model (they
        # were 0, Meff2 included, except for RPH and Horndeski, whose own values the same formulas reproduce). Pure EFT
        # Omega = a: Meff2 = 1 + a, alphaM = a/(1 + a), alphaB = a/(2 (1 + a)) (EFTCAMB convention), alphaT = 0
        a = np.array([0.1, 0.5, 0.9])
        cp = self.read("2_PEFT_Omega_lin_1")
        res = camb.get_results(cp)
        for f in (cp.EFTCAMB.get_eft_functions(res, a)[1], cp.EFTCAMB.get_scale_evolution(res, [1e-2], a)[1][0]):
            af = np.ravel(f["a"])  # the evolution is sampled in conformal time: a to ~1e-10
            np.testing.assert_allclose(np.ravel(f["Meff2"]), 1 + af, rtol=1e-10)
            np.testing.assert_allclose(np.ravel(f["alphaM"]), af / (1 + af), rtol=1e-10)
            np.testing.assert_allclose(np.ravel(f["alphaB"]), af / (2 * (1 + af)), rtol=1e-10)
            np.testing.assert_allclose(np.ravel(f["alphaT"]), 0, atol=1e-15)

    def test_mass_eigenvalues(self):
        # compute_mass_factors (arXiv:1609.03599 App. C, corrected: d^2 EFTc/da^2, the matter rho + p as dust with a
        # consistent Hdotdot, quadruple precision): EFT_mu1 is the scalaron mass ((1 + f_R)/f_RR - R)/3 of a heavy
        # f(R) scalaron (it was +-1e13 H^2 or inf) and of order H^2 for quintessence (it was +-1e15 H^2)
        a = np.array([0.5, 1.0])
        mu1, scalaron = self.fR_mass("EFTB0 = 1e-4", a)
        np.testing.assert_allclose(mu1, scalaron, rtol=1e-2)
        cp = self.read("5_quint_exponential")
        ev = cp.EFTCAMB.get_scale_evolution(camb.get_results(cp), [1e-3], a)[1][0]
        self.assertTrue(np.all(np.abs(np.ravel(ev["EFT_mu1"])) < 5 * (np.ravel(ev["adotoa"]) / a) ** 2))

    def test_mass_stability(self):
        # EFT_mass_stability rejects a tachyonic f(R) scalaron (B0 < 0 with the f_RR prior off: mu_1 ~ -4e7 H^2); it
        # was experimental and ignored, and the run then failed in the integrator. Healthy models still run.
        with self.assertRaisesRegex(CAMBError, "theory unstable"):
            camb.get_results(
                self.read("4_Des_F_R_LCDM_1", "EFTB0 = -1e-4\nEFT_additional_priors = F\nEFT_mass_stability = T")
            )
        for model, extra in (("4_Des_F_R_LCDM_1", "EFTB0 = 1e-4"), ("5_quint_exponential", "")):
            with self.subTest(model=model):
                cp = self.read(model, extra + "\nEFT_mass_stability = T")
                self.assertGreater(camb.get_results(cp).get_sigma8_0(), 0)

    def test_mass_stability_methods(self):
        # EFT_mass_stability_method = 2: the normal-form mass of pi (the f(R) scalaron mass, also where method 1 is not
        # resolved), tested only where pi itself also grows faster than the rate. Hilltop quintessence (axion potential,
        # n = 1, m = 8: V'' down to -58 H^2) is refused by method 1 at rate 3, but pi = delta phi/phi' grows no faster
        # than with m = 1 (delta phi follows the rolling background): method 2 accepts it. Both refuse a tachyonic f(R).
        a = np.array([0.1, 0.5, 1.0])
        mu1, scalaron = self.fR_mass("EFTB0 = 1e-2\nEFT_mass_stability_method = 2", a)
        np.testing.assert_allclose(mu1, scalaron, rtol=1e-4)
        hilltop = "n = 1\nm = 8\nEFT_mass_stability = T\nEFT_mass_stability_rate = 3\nEFT_mass_stability_method = "
        with self.assertRaisesRegex(CAMBError, "theory unstable"):
            camb.get_results(self.read("5_quint_axion", hilltop + "1"))
        self.assertGreater(camb.get_results(self.read("5_quint_axion", hilltop + "2")).get_sigma8_0(), 0)
        tachyon = "EFTB0 = -1e-4\nEFT_additional_priors = F\nEFT_mass_stability = T\nEFT_mass_stability_method = 2"
        with self.assertRaisesRegex(CAMBError, "theory unstable"):
            camb.get_results(self.read("4_Des_F_R_LCDM_1", tachyon))
        with self.assertRaisesRegex(CAMBValueError, "invalid EFT_mass_stability_method"):
            self.read("5_quint_axion", "EFT_mass_stability_method = 3")

    def test_3p1_b22(self):
        # 3+1 designer gradient test with the corrected B22 (+ b3/2, missing in arXiv:2508.01378 App. D): this healthy
        # model was refused by the model's own gradient test
        model = (
            "EFTflag = 3\nDesignerEFTmodel = 3\nEFTwDE = 1\nEFTw0 = -0.95\nDesigner3p1LuminalGW = T\n"
            "Designer3p1modelsigma1 = 1\nDesigner3p1sigma10 = 1.0\nDesigner3p1modelsigma1_ODE = 2\n"
            "Designer3p1sigma1_ODE0 = -0.5\nDesigner3p1modelrho1N = 0\nDesigner3p1modelrho2 = 1\nDesigner3p1rho20 = 1.0\n"
            "EFT_ghost_stability = F\nEFT_gradient_stability = F\nEFT_additional_priors = T"
        )
        self.assertGreater(camb.get_results(self.read("2_PEFT_wDE_wCDM_1", model)).get_sigma8_0(), 0)

    def test_3p1_dust_ghost(self):
        # with CDM the 3+1 scalar is a ghost at every k for sigma1 > 1 (alpha_L < 0): the model's own test used the
        # vacuum conditions and accepted this one (sigma8 = 5.7e8); sigma1 < 1 still runs
        model = (
            "EFTflag = 3\nDesignerEFTmodel = 3\nEFTwDE = 1\nEFTw0 = -0.95\nDesigner3p1LuminalGW = T\n"
            "Designer3p1modelsigma1 = 1\nDesigner3p1sigma10 = 1.0\nDesigner3p1modelsigma1_ODE = 2\n"
            "Designer3p1modelrho1N = 0\nDesigner3p1modelrho2 = 1\nDesigner3p1rho20 = 1.e-3\n"
            "EFT_ghost_stability = F\nEFT_gradient_stability = F\nEFT_additional_priors = T\n"
        )
        with self.assertRaisesRegex(CAMBError, "theory unstable"):
            camb.get_results(self.read("2_PEFT_wDE_wCDM_1", model + "Designer3p1sigma1_ODE0 = 1.e-7"))
        res = camb.get_results(self.read("2_PEFT_wDE_wCDM_1", model + "Designer3p1sigma1_ODE0 = -1.e-7"))
        self.assertGreater(res.get_sigma8_0(), 0)

    def test_horndeski_coefficient_identities(self):
        # the generated Horndeski background coefficients are conformal-time derivatives of each other: checked at
        # random inputs with central differences along a' = a^2 H, X = dphi^2/(2 a^2), dG/dtau = G_phi phi' + G_X X'
        # (a G4pxx for G4ppp typo broke cs1p, d0 and dp0)
        from camb.baseconfig import camblib

        g_names = [
            "G2",
            "G2p",
            "G2x",
            "G2pp",
            "G2px",
            "G2xx",
            "G2ppp",
            "G2ppx",
            "G2pxx",
            "G2xxx",
            "G3p",
            "G3x",
            "G3pp",
            "G3px",
            "G3xx",
            "G3ppp",
            "G3ppx",
            "G3pxx",
            "G3xxx",
            "G4",
            "G4p",
            "G4x",
            "G4pp",
            "G4px",
            "G4xx",
            "G4ppp",
            "G4ppx",
            "G4pxx",
            "G4xxx",
            "G4pppp",
            "G4pppx",
            "G4ppxx",
            "G4pxxx",
            "G4xxxx",
            "G5",
            "G5p",
            "G5x",
            "G5pp",
            "G5px",
            "G5xx",
            "G5ppp",
            "G5ppx",
            "G5pxx",
            "G5xxx",
            "G5pppp",
            "G5pppx",
            "G5ppxx",
            "G5pxxx",
            "G5xxxx",
        ]
        signatures = {
            "backeq": ("horndeski_backeq_coefficients", ("a", "H", "dphi", "rptot"), 6),
            "backeq_p": ("horndeski_backeq_coefficients_p", ("a", "H", "dH", "dphi", "ddphi", "rptot", "drptot"), 6),
            "designer": ("horndeski_designereq_coefficients", ("a", "H", "dH", "ddH", "dphi", "rptot", "drptot"), 3),
            "designer_p": (
                "horndeski_designereq_coefficients_p",
                ("a", "H", "dH", "ddH", "dddH", "dphi", "ddphi", "rptot", "drptot", "ddrptot"),
                2,
            ),
        }

        def call(name, v):
            symbol, first, n_out = signatures[name]
            outputs = [c_double() for _ in range(n_out)]
            getattr(camblib, symbol)(*[c_double(v[k]) for k in first + tuple(g_names)], *[byref(o) for o in outputs])
            return np.array([o.value for o in outputs])

        def rates(v):
            xdot = v["dphi"] * v["ddphi"] / v["a"] ** 2 - v["dphi"] ** 2 * v["H"] / v["a"]
            r = {"a": v["a"] ** 2 * v["H"], "H": v["dH"], "dH": v["ddH"], "ddH": v["dddH"], "dphi": v["ddphi"],
                 "ddphi": v["dddphi"], "rptot": v["drptot"], "drptot": v["ddrptot"]}  # fmt: skip
            for g in g_names:
                base, n_p, n_x = g[:2], g[2:].count("p"), g[2:].count("x")
                g_p, g_x = base + "p" * (n_p + 1) + "x" * n_x, base + "p" * n_p + "x" * (n_x + 1)
                r[g] = v.get(g_p, 0.0) * v["dphi"] + v.get(g_x, 0.0) * xdot
            return r

        def shifted(v, s, h=1e-6):
            r = rates(v)
            return {k: x + s * h * r.get(k, 0.0) for k, x in v.items()}

        rng = np.random.default_rng(1)
        for _ in range(10):
            v = {k: rng.uniform(-1, 1) for k in g_names}
            v.update(a=rng.uniform(0.3, 1), H=rng.uniform(0.5, 2), dphi=rng.uniform(0.2, 1), rptot=rng.uniform(0, 1))
            v.update({k: rng.uniform(-1, 1) for k in ("dH", "ddH", "dddH", "ddphi", "dddphi", "drptot", "ddrptot")})
            # backeq_p = d/dtau backeq
            fd = (call("backeq", shifted(v, 1)) - call("backeq", shifted(v, -1))) / 2e-6
            np.testing.assert_allclose(call("backeq_p", v), fd, rtol=1e-6)
            # d0 + d1 phi'' + d2 phi''^2 = d/dtau(cp1 phi'' + ch1 H' + cs1) - cp1 phi'''
            d, b, bp = call("designer", v), call("backeq", v), call("backeq_p", v)
            lhs = d[0] + d[1] * v["ddphi"] + d[2] * v["ddphi"] ** 2
            self.assertAlmostEqual(
                lhs / (bp[0] * v["ddphi"] + bp[1] * v["dH"] + b[1] * v["ddH"] + bp[2]), 1, delta=1e-10
            )
            # dp0 + dp1 phi''' = d/dtau(d0 + d1 phi'' + d2 phi''^2) where the designer equation holds, cp1 = cp1' = 0
            # (G4p and G4pp enter cp1 and cp1' linearly)
            c0, c1 = (call("backeq", {**v, "G4p": x})[0] for x in (0.0, 1.0))
            v["G4p"] = -c0 / (c1 - c0)
            c0, c1 = (call("backeq_p", {**v, "G4pp": x})[0] for x in (0.0, 1.0))
            v["G4pp"] = -c0 / (c1 - c0)

            def poly(w):
                dd = call("designer", w)
                return dd[0] + dd[1] * w["ddphi"] + dd[2] * w["ddphi"] ** 2

            dp = call("designer_p", v)
            self.assertAlmostEqual(
                (dp[0] + dp[1] * v["dddphi"]) / ((poly(shifted(v, 1)) - poly(shifted(v, -1))) / 2e-6), 1, delta=1e-6
            )

    def test_horndeski_pi_path_turning_point(self):
        # on the pi path (EFTCAMB_evolve_delta_phi = F) pi = delta_phi/phi' is singular where the Horndeski field velocity
        # changes sign: 5_hdsk_TG_3 (a = 1.69e-3) ran with sigma8 1.77 instead of 0.903 and zero CMB C_l, even with the
        # stability check on. It is refused whatever the stability flags; TG_2 (pi path, no crossing) still runs.
        message = r"field velocity changes sign at a = 1\.692E-03"
        for extra in ("", "\nEFTCAMB_skip_stability = F"):
            with self.subTest(extra=extra), self.assertRaisesRegex(CAMBError, message):
                camb.get_background(self.read("5_hdsk_TG_3", "EFTCAMB_evolve_delta_phi = F" + extra))
        self.assertGreater(camb.get_background(self.read("5_hdsk_TG_2")).hubble_parameter(0), 0)

    def test_horndeski_gr_limit(self):
        # LCDM written as Horndeski (G2 = const, G4 = 1/2): rho_DE + p_DE comes from the background equation without
        # subtracting matter, so c and gamma1 are exactly 0 (they were noise up to 1e-4 and the pi path aborted)
        model = "Horndeski_freefunc3_model = 0\nHorndeskic00 = 1.1653675624178562e-07\n"
        model += "Horndeskic01 = 0\nEFTCAMB_skip_RGR = F"
        a = np.logspace(-4, 0, 9)
        sigma8 = {}
        for dphi in ("T", "F"):
            cp = self.read("5_hdsk_TG_1", f"{model}\nEFTCAMB_evolve_delta_phi = {dphi}")
            res = camb.get_results(cp)
            ev = cp.EFTCAMB.get_scale_evolution(res, [0.05], a)[1][0]
            self.assertFalse(np.any(ev["EFTc"]) or np.any(ev["EFTGamma1V"]))
            sigma8[dphi] = res.get_sigma8_0()
        self.assertAlmostEqual(sigma8["F"] / sigma8["T"], 1, delta=1e-12)

    def test_horava_domain(self):
        # Horava (code variables lambda = lambda_p - 1, xi = xi_p - 1): 0 < eta < 2 xi + 2 and lambda >= 0, except at
        # the GR point (arXiv:1508.01787 Sec. 2.4)
        for label, extra in (
            ("eta = 0", "Horava_eta = 0"),
            ("eta = 2 + 2 xi", "Horava_eta = 1.998"),
            ("lambda = -2/3", "Horava_lambda = -0.6666666666666666"),
            ("lambda = -1", "Horava_lambda = -1"),
            ("lambda = -0.1", "Horava_lambda = -0.1"),
        ):
            with self.subTest(refused=label), self.assertRaises(CAMBError):
                camb.get_results(self.read("5_Horava_eta_3", extra))
        for label, extra in (
            ("GR point", "Horava_eta = 0\nHorava_lambda = 0\nHorava_xi = 0"),
            ("lambda = 0", "Horava_lambda = 0"),
        ):
            with self.subTest(accepted=label):
                self.assertGreater(camb.get_results(self.read("5_Horava_eta_3", extra)).get_sigma8_0(), 0)

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

    def test_beyond_today_refused(self):
        # EFT results stop today (a = 1): requests with z < 0, a > 1 or eta > tau0 are refused
        res = camb.get_results(self.read("2_PEFT_wDE_wCDM_1"))
        eft = res.Params.EFTCAMB
        requests = {
            "hubble_parameter": lambda: res.hubble_parameter(-0.1),
            "angular_diameter_distance": lambda: res.angular_diameter_distance([0.5, -0.1]),
            "angular_diameter_distance2": lambda: res.angular_diameter_distance2(-0.1, 1.0),
            "comoving_radial_distance": lambda: res.comoving_radial_distance(np.array([-0.1, 1.0])),
            "luminosity_distance": lambda: res.luminosity_distance(-0.1),
            "conformal_time": lambda: res.conformal_time(-0.1),
            "sound_horizon": lambda: res.sound_horizon(-0.1),
            "physical_time": lambda: res.physical_time(-0.1),
            "conformal_time_a1_a2": lambda: res.conformal_time_a1_a2(0.5, 1.2),
            "get_background_densities": lambda: res.get_background_densities(1.5),
            "get_dark_energy_rho_w": lambda: res.get_dark_energy_rho_w(1.5),
            "get_Omega": lambda: res.get_Omega("de", -0.1),
            "get_background_redshift_evolution": lambda: res.get_background_redshift_evolution([-0.1]),
            "get_background_time_evolution": lambda: res.get_background_time_evolution([res.tau0 * 1.01]),
            "redshift_at_comoving_radial_distance": lambda: res.redshift_at_comoving_radial_distance(-10.0),
            "get_redshift_evolution": lambda: res.get_redshift_evolution(0.1, [-0.1]),
            "get_eft_functions": lambda: eft.get_eft_functions(res, [0.5, 1.2]),
            "get_evolution": lambda: eft.get_evolution(res, 0.1, [100.0, res.tau0 * 1.01]),
        }
        for name, request in requests.items():
            with self.subTest(request=name), self.assertRaises(CAMBValueError):
                request()
        # today is allowed
        self.assertAlmostEqual(res.hubble_parameter(0.0) / res.Params.H0, 1, delta=1e-6)
        self.assertAlmostEqual(res.get_dark_energy_rho_w(1.0)[0], 1, delta=1e-6)
        self.assertEqual(eft.get_eft_functions(res, [0.5, 1.0])[1].shape, (2,))
        # z_outputs are refused by the Fortran parameter check
        with self.assertRaises(CAMBError):
            res.get_BAO([0.5, -0.1], res.Params)
        # plain CAMB is not restricted
        self.assertGreater(camb.get_background(camb.set_params(H0=67.0)).hubble_parameter(-0.1), 0)
        # a = 1 through conformal times: conformal_time(z = 0) of RPH with 6000-point tables (the former default)
        # exceeds tau0 by ~2e-8 for this array
        cp = self.read("3_RPH_KT_lin_1", "model_background_num_points = 6000")
        res = camb.get_results(cp)
        a = np.array([0.05, 0.1, 0.3, 0.6, 1.0])
        self.assertGreater(res.conformal_time(1 / a - 1)[-1], res.tau0)
        self.assertEqual(cp.EFTCAMB.get_scale_evolution(res, [0.05, 0.5], a)[1].shape, (2, 5))
        self.assertEqual(res.get_background_redshift_evolution(1 / a - 1, format="array").shape[0], 5)

    def test_wDE_prior(self):
        # the phenomenological prior w_DE <= -1/3 is off by default and refuses with its own message when selected
        for model in ("2_PEFT_wDE_wCDM_1", "4_Des_mc5e_wCDM_1"):
            with self.subTest(model=model):
                cp = self.read(model, "EFTw0 = -0.3")
                self.assertFalse(cp.EFTCAMB.EFT_wDE_prior)
                self.assertGreater(camb.get_results(cp).get_sigma8_0(), 0)
                cp = self.read(model, "EFTw0 = -0.3\nEFT_wDE_prior = T")
                # the fields next to the new flag keep their values (ctypes layout)
                self.assertTrue(cp.EFTCAMB.EFT_wDE_prior)
                self.assertTrue(cp.EFTCAMB.EFT_additional_priors)
                self.assertFalse(cp.EFTCAMB.EFT_positivity_bounds)
                with self.assertRaisesRegex(CAMBError, "w_DE > -1/3 at a = .*EFT_wDE_prior"):
                    camb.get_results(cp)
                # the ini round trip keeps the flag
                name = os.path.join(self.tmp, f"written_{model}.ini")
                cp.write_ini(name)
                self.assertTrue(camb.read_ini(name).EFTCAMB.EFT_wDE_prior)
                # set_params selects it too
                eft = {**cp.EFTCAMB.read_parameters(), "EFTw0": -0.3}
                with self.assertRaisesRegex(CAMBError, "EFT_wDE_prior"):
                    camb.get_results(camb.set_params(H0=67.0, lmax=600, **eft))
                # an accelerating w_DE passes the prior
                cp = self.read(model, "EFTw0 = -0.9\nEFT_wDE_prior = T")
                self.assertGreater(camb.get_results(cp).get_sigma8_0(), 0)

    def test_no_kinetic_term_refused(self):
        # a pi field with no kinetic term (A = 0) and a sourced equation is refused, not set to zero
        no_kinetic = "PureEFTmodelGamma1 = 0"
        with self.assertRaisesRegex(CAMBError, r"no kinetic term \(A = 0\) at a ="):
            camb.get_results(self.read("2_PEFT_gamma3_lin_1", no_kinetic))
        # without the stability check the evolution refuses it
        with self.assertRaisesRegex(CAMBError, r"no kinetic term \(A = 0\) but a sourced equation"):
            camb.get_results(self.read("2_PEFT_gamma3_lin_1", no_kinetic + "\nEFTCAMB_skip_stability = T"))
        # with a kinetic term the model runs
        self.assertGreater(camb.get_results(self.read("2_PEFT_gamma3_lin_1")).get_sigma8_0(), 0)
        # GR has no pi kinetic term either, and an empty pi equation: pi = 0 (here evolved, skipping the return to GR)
        sigma8 = camb.get_results(self.read("2_PEFT_wDE_LCDM_1")).get_sigma8_0()
        sigma8_evolved = camb.get_results(self.read("2_PEFT_wDE_LCDM_1", "EFTCAMB_skip_RGR = T")).get_sigma8_0()
        self.assertAlmostEqual(sigma8_evolved / sigma8, 1, delta=1e-4)

    def test_rph_alphaM_feedback(self):
        # the feedback of an RPH alpha_M(a) model prints the real RPH_M0 (an integer format aborted the process):
        # run it in a subprocess, so that an abort fails the test instead of the test runner
        # (at the default feedback level: the fast settings switch the feedback off)
        ini = os.path.join(self.tmp, "rph_alphaM.ini")
        with open(ini, "w") as f:
            f.write(f"DEFAULT({os.path.join(PARAMS, '3_RPH_MB_OmegaDE_1.ini')})\n")
            f.write("RPHalphaMmodel = 1\nRPHalphaM0 = 0.01\nRPHalphaMmodel_ODE = 0\n")
        code = (
            f"import camb; cp = camb.read_ini({ini!r}); cp.EFTCAMB.feedback(print_params=True); print('process alive')"
        )
        r = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
            env=dict(os.environ, PYTHONPATH=ROOT, OMP_NUM_THREADS="1"),
        )
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        self.assertIn("RPH_M0", r.stdout)
        self.assertIn("process alive", r.stdout)

    def test_rph_gr_point(self):
        # RPH reads H, its derivatives and rho_DE from linearly interpolated tables: with 6000 points the GR point (all
        # alpha = 0, EFT never on) was off the Pure EFT LCDM by TT 9.8e-6, sigma8 1.9e-6; the 60000-point default 1e-7
        tt = lambda res: res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][2:, 0]
        lcdm = camb.get_results(self.read("2_PEFT_wDE_LCDM_1"))
        res = camb.get_results(self.read("3_RPH_KT_lin_1", "RPHkineticitymodel = 0\nRPHtensormodel = 0"))
        self.assertGreater(res.Params.EFTCAMB.EFTCAMB_pert_turn_on, 1)
        np.testing.assert_allclose(tt(res), tt(lcdm), rtol=1e-6)
        self.assertAlmostEqual(res.get_sigma8_0() / lcdm.get_sigma8_0(), 1, delta=2e-7)

    def test_rph_wde_normalisation(self):
        # RPH with w_DE != -1: rho_DE(a = 1) = 3 H0^2 Omega_DE, so H(z = 0) = H0, with both integration branches
        z = np.array([0.0, 0.5, 2.0, 10.0])
        for wde in ("RPHwDE = 1\nRPHw0 = -0.9", "RPHwDE = 2\nRPHw0 = -0.9\nRPHwa = 0.1"):
            with self.subTest(wde=wde):
                cp = self.read("3_RPH_MK_const_1", wde)
                h = camb.get_background(cp).hubble_parameter(z)
                self.assertAlmostEqual(h[0] / cp.H0, 1, delta=1e-9)
                h_today = camb.get_background(
                    self.read("3_RPH_MK_const_1", wde + "\nRPHintegratefromtoday = T")
                ).hubble_parameter(z)
                np.testing.assert_allclose(h, h_today, rtol=1e-9)

    def theta_ini(self, model, theta_line):
        # the model ini with hubble removed from its DEFAULT chain and a theta input instead
        def copy(path):
            lines = []
            with open(path) as f:
                for line in f:
                    if line.strip().startswith("DEFAULT("):
                        lines.append(f"DEFAULT({copy(os.path.join(PARAMS, line.strip()[8:-1]))})\n")
                    elif line.split("=")[0].strip() != "hubble":
                        lines.append(line)
            name = os.path.join(self.tmp, "nohubble_" + os.path.basename(path))
            with open(name, "w") as f:
                f.writelines(lines)
            return name

        name = os.path.join(self.tmp, f"{model}_theta_{abs(hash(theta_line))}.ini")
        with open(name, "w") as f:
            f.write(f"DEFAULT({copy(os.path.join(PARAMS, model + '.ini'))})\n{FAST}output_root = {self.tmp}/out\n")
            f.write(f"EFTCAMB_feedback_level = 0\n{theta_line}\n")
        return name

    def test_theta_input(self):
        # Horndeski without designer potential sets H0 from its own background (also with shooting): theta refused
        for extra in ("", "Horndeski_shooting = T"):
            with self.subTest(shooting=extra):
                for key in ("thetastar", "cosmomc_theta"):
                    with self.assertRaisesRegex(CAMBValueError, "H0 .* is set by the EFT background"):
                        camb.read_ini(self.theta_ini("5_hdsk_TG_1", f"{extra}\n{key} = 0.0104"))
                cp = self.read("5_hdsk_TG_1", f"EFTCAMB_feedback_level = 0\n{extra}")
                with self.assertRaisesRegex(CAMBValueError, "set by the EFT background"):
                    cp.set_cosmology(H0=None, thetastar=0.0104, ombh2=cp.ombh2, omch2=cp.omch2)
                with self.assertRaisesRegex(CAMBValueError, "set by the EFT background"):
                    cp.set_H0_for_theta(0.0104, cosmomc_approx=True)
        # other EFT models keep the theta -> H0 solve on both paths
        cp = camb.read_ini(self.theta_ini("4_Des_F_R_LCDM_1", "thetastar = 0.0104"))
        self.assertAlmostEqual(camb.get_background(cp).get_derived_params()["thetastar"] / 1.04, 1, delta=1e-5)
        cp = self.read("4_Des_F_R_LCDM_1")
        cp.set_H0_for_theta(0.0104)
        self.assertAlmostEqual(camb.get_background(cp).get_derived_params()["thetastar"] / 1.04, 1, delta=1e-5)

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

    def test_test_inis_keys_read(self):
        # every EFTCAMB key a test ini sets is read by its model: a typo (EFtw2) or a key of another parametrization
        # (an *exp key with a linear model, EFTGamma40 in the gamma3 Horndeski mode) is silently ignored by the reader
        unread = {}
        for path in sorted(glob.glob(os.path.join(PARAMS, "*.ini"))):
            with open(path) as f:
                text = f.read()
            if not re.search(r"^\s*EFTflag\s*=\s*[1-9]", text, re.MULTILINE):
                continue
            keys = re.findall(r"^\s*([A-Za-z_]\w*)\s*=", text, re.MULTILINE)
            model = os.path.basename(path)[:-4]
            cp = self.read(model)  # keep the parent of cp.EFTCAMB alive
            read = set(cp.EFTCAMB.read_parameters())
            missing = [k for k in keys if eftcamb.is_eftcamb_parameter(k) and k not in read]
            missing = [k for k in missing if k not in ("output_root", "feedback_level")]
            if missing:
                unread[model] = missing
        self.assertEqual(unread, {})

    def test_regression(self):
        for model, (sigma8, tt) in REFERENCE.items():
            with self.subTest(model=model):
                res = camb.get_results(self.read(model))
                cl = res.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][[10, 100, 400], 0]
                self.assertAlmostEqual(res.get_sigma8_0() / sigma8, 1, delta=1e-5)
                np.testing.assert_allclose(cl, tt, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
