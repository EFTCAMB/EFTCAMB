"""Worker for track 8c checks; runs inside a subprocess started by t8c_lib.run_task.

cfg = {tree, ini, task, args, out}.  Tasks:
  generic : read_ini, then optional Python-side setup from args, then outputs listed in args["outputs"].
  setparams : build CAMBparams with camb.set_params(**args["params"]) (pure Python path, no ini).
"""

from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np


def _set_attr(obj, path: str, value) -> None:
    parts = path.split(".")
    for p in parts[:-1]:
        obj = getattr(obj, p)
    setattr(obj, parts[-1], value)


Z_BG = np.concatenate([np.linspace(0, 10, 101), np.logspace(1.05, 4, 40)])


def setup(camb, cp, args, meta):
    from camb import model

    for name, value in args.get("config", {}).items():
        # global camb.config settings (e.g. AccuracyTarget, read only by the command-line driver from an ini)
        setattr(camb.config, name, value)
        meta[f"config_{name}"] = getattr(camb.config, name)
    for path, value in args.get("set_attrs_pre", {}).items():
        _set_attr(cp, path, value)
    if "dark_energy" in args:
        de = dict(args["dark_energy"])
        name = de.pop("model", None)
        if name:
            cp.set_classes(dark_energy_model=name)
        cp.DarkEnergy.set_params(**de)
    if "cosmology" in args:
        cp.set_cosmology(**args["cosmology"])
    if "set_for_lmax" in args:
        cp.set_for_lmax(**args["set_for_lmax"])
    if "set_H0_for_theta" in args:
        cp.set_H0_for_theta(**args["set_H0_for_theta"])
    if "theta_via_set_cosmology" in args:
        # re-set the current cosmology through set_cosmology (EFTCAMB-modified H0 search range) with theta
        kind, value = args["theta_via_set_cosmology"]
        omnu = float(cp.omnuh2)
        nmass = int(cp.num_nu_massive)
        cp.set_cosmology(H0=None, ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, mnu=None if omnu > 0 else 0.0,
                         omnuh2_active=omnu if omnu > 0 else None, num_massive_neutrinos=max(nmass, 1),
                         nnu=cp.N_eff, YHe=cp.YHe, TCMB=cp.TCMB, **{kind: value})
    if "set_classes_nl" in args:
        cp.set_classes(non_linear_model=args["set_classes_nl"])
    if "NonLinear" in args:
        cp.NonLinear = getattr(model, args["NonLinear"])
    if "nonlinear_model" in args:
        cp.NonLinearModel.set_params(**args["nonlinear_model"])
    if "matter_power" in args:
        cp.set_matter_power(**args["matter_power"])
    if "source_windows" in args:
        from camb.sources import GaussianSourceWindow

        cp.SourceWindows = [GaussianSourceWindow(**w) for w in args["source_windows"]]
    if "source_terms" in args:
        for k, v in args["source_terms"].items():
            setattr(cp.SourceTerms, k, v)
    for path, value in args.get("set_attrs", {}).items():
        _set_attr(cp, path, value)
    meta["NonLinear"] = str(cp.NonLinear)
    meta["max_l"] = int(cp.max_l)
    meta["max_eta_k"] = float(cp.max_eta_k)
    meta["kmax"] = float(cp.Transfer.kmax)
    meta["DarkEnergy"] = type(cp.DarkEnergy).__name__
    meta["H0"] = float(cp.H0)
    meta["omk"] = float(cp.omk)
    meta["EFTflag"] = int(cp.EFTCAMB.EFTflag) if hasattr(cp, "EFTCAMB") else None


def outputs(camb, cp, args, out, meta):
    outs = args.get("outputs", ["cmb"])
    t0 = time.time()
    if "transfer_only" in outs:
        results = camb.get_transfer_functions(cp)
        if "power_from_transfer" in outs:
            results.power_spectra_from_transfer()
    elif "background_only" in outs:
        results = camb.get_background(cp)
    else:
        results = camb.get_results(cp)
    meta["time_results"] = time.time() - t0
    meta["H0_final"] = float(results.Params.H0)
    try:
        meta["derived"] = {k: float(v) for k, v in results.get_derived_params().items()}
    except Exception as e:  # noqa: BLE001
        meta["errors"]["derived"] = repr(e)
    if "cmb" in outs:
        lmax = args.get("cl_lmax")
        cls = results.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True, lmax=lmax)
        for key in ("total", "unlensed_scalar", "lensed_scalar", "lens_potential", "tensor"):
            if key in cls:
                out[f"cl_{key}"] = np.asarray(cls[key])
    if "sources" in outs:
        d = results.get_source_cls_dict(raw_cl=True)
        keys = sorted(d)
        out["src_keys"] = np.array(keys)
        out["src_cls"] = np.array([d[k] for k in keys])
    if "pk_lin" in outs or "pk_nl" in outs:
        k, z, pk = results.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False)
        out["pk_k"], out["pk_z"], out["pk"] = np.asarray(k), np.asarray(z), np.asarray(pk)
        out["sigma8"] = np.atleast_1d(results.get_sigma8())
    if "pk_nl" in outs:
        k, z, pk = results.get_nonlinear_matter_power_spectrum(hubble_units=False, k_hunit=False)
        out["pknl_k"], out["pknl"] = np.asarray(k), np.asarray(pk)
    if "pk_weyl" in outs:
        k, z, pk = results.get_linear_matter_power_spectrum(var1="Weyl", var2="Weyl", hubble_units=False,
                                                            k_hunit=False)
        out["pkweyl"] = np.asarray(pk)
        k, z, pk = results.get_linear_matter_power_spectrum(var1="delta_tot", var2="Weyl", hubble_units=False,
                                                            k_hunit=False)
        out["pkdw"] = np.asarray(pk)
    if "transfer" in outs:
        td = results.get_matter_transfer_data()
        out["tr_q"] = np.asarray(td.q)
        out["tr"] = np.asarray(td.transfer_data)
        out["tr_z"] = np.asarray(td.transfer_z, dtype=float)
    if "background" in outs:
        a = 1 / (1 + Z_BG)
        out["bg_z"] = Z_BG
        out["bg_H"] = np.asarray(results.hubble_parameter(Z_BG))
        out["bg_DA"] = np.asarray(results.angular_diameter_distance(Z_BG))
        out["bg_tau"] = np.asarray(results.conformal_time(Z_BG))
        dens = results.get_background_densities(a, format="array")
        out["bg_dens"] = np.asarray(dens)
        rho, w = results.get_dark_energy_rho_w(a)
        out["bg_de_rho"], out["bg_de_w"] = np.asarray(rho), np.asarray(w)
        out["bg_Omega_de"] = np.asarray(results.get_Omega("de", Z_BG))
        out["bg_physical_time"] = np.asarray(results.physical_time(Z_BG[Z_BG < 1000]))
        try:
            ev = results.get_background_redshift_evolution(Z_BG[Z_BG < 3000], format="array")
            out["bg_evol"] = np.asarray(ev)
        except Exception as e:  # noqa: BLE001
            meta["errors"]["bg_evol"] = repr(e)
        meta["age_results"] = float(results.physical_time(0.0))
        meta["age_get_age"] = float(camb.get_age(cp))
        meta["tau0"] = float(results.tau0)
        meta["cosmomc_theta"] = float(results.cosmomc_theta())
        if hasattr(cp, "EFTCAMB") and cp.EFTCAMB.EFTflag != 0:
            try:
                fields, arr = cp.EFTCAMB.get_eft_functions(results, a[Z_BG < 1000])
                for f in ("adotoa", "grhov_t", "grhoq", "gpresq", "grhom_t", "EFTc", "EFTLambda", "EFTOmegaV"):
                    out[f"eft_{f}"] = np.ascontiguousarray(arr[f])
            except Exception as e:  # noqa: BLE001
                meta["errors"]["eft_functions"] = repr(e)
    if "zre" in outs:
        for tau in args.get("zre_taus", [0.054, 0.09]):
            meta[f"zre_from_tau_{tau}"] = float(camb.get_zre_from_tau(cp, tau))
        meta["reion_redshift"] = float(results.Params.Reion.redshift)
        meta["reion_tau"] = float(results.Params.Reion.optical_depth)
        meta["reion_tau_derived_zre"] = float(results.get_derived_params().get("zrei", np.nan))
    if "theta" in outs:
        meta["cosmomc_theta"] = float(results.cosmomc_theta())
        meta["thetastar"] = float(results.get_derived_params()["thetastar"]) / 100
    if "weyl_limber" in outs:
        # CAMBdemo recipe: Limber C_L^phiphi from P_Weyl(k,z) vs CAMB's line-of-sight lensing potential
        from camb import model as cmodel

        lmax = 2000
        nz = 200
        chistar = results.conformal_time(0) - results.tau_maxvis
        chis = np.linspace(0, chistar, nz)
        zs = results.redshift_at_comoving_radial_distance(chis)
        dchis = (chis[2:] - chis[:-2]) / 2
        chis, zs = chis[1:-1], zs[1:-1]
        kmax = 10.0
        PK = camb.get_matter_power_interpolator(cp, nonlinear=False, hubble_units=False, k_hunit=False, kmax=kmax,
                                                var1=cmodel.Transfer_Weyl, var2=cmodel.Transfer_Weyl, zmax=zs[-1])
        win = ((chistar - chis) / (chis**2 * chistar)) ** 2
        ls = np.arange(2, lmax + 1, dtype=np.float64)
        cl_kappa = np.zeros(ls.shape)
        w = np.ones(chis.shape)
        for i, ell in enumerate(ls):
            kk = (ell + 0.5) / chis
            w[:] = 1
            w[kk < 1e-4] = 0
            w[kk >= kmax] = 0
            cl_kappa[i] = np.dot(dchis, w * PK.P(zs, kk, grid=False) * win / kk**4)
        cl_kappa *= (ls * (ls + 1)) ** 2
        out["limber_ls"] = ls
        out["limber_pp"] = 4 * cl_kappa / 2 / np.pi
        out["los_pp"] = np.asarray(results.get_lens_potential_cls(lmax)[2:, 0])
    if "nu" in outs:
        P = results.Params
        nu = {"omnuh2": float(P.omnuh2), "num_nu_massive": int(P.num_nu_massive),
              "nu_mass_eigenstates": int(P.nu_mass_eigenstates),
              "nu_mass_degeneracies": [float(x) for x in P.nu_mass_degeneracies],
              "nu_mass_fractions": [float(x) for x in P.nu_mass_fractions],
              "nu_masses": [float(x) for x in results.nu_masses], "grhormass": [float(x) for x in results.grhormass]}
        if hasattr(P, "EFTCAMB") and P.EFTCAMB.EFTflag != 0:
            try:
                c = P.EFTCAMB_parameter_cache
                nu["cache_nu_masses"] = [float(x) for x in c.nu_masses]
                nu["cache_grhormass"] = [float(x) for x in c.grhormass]
                nu["cache_Nu_mass_eigenstates"] = int(c.Nu_mass_eigenstates)
                n = nu["nu_mass_eigenstates"]
                nu["cache_match"] = bool(np.allclose(nu["cache_nu_masses"][:n], nu["nu_masses"][:n], rtol=1e-14)
                                         and np.allclose(nu["cache_grhormass"][:n], nu["grhormass"][:n], rtol=1e-14)
                                         and nu["cache_Nu_mass_eigenstates"] == n)
            except Exception as e:  # noqa: BLE001
                meta["errors"]["nu_cache"] = repr(e)
        meta["nu"] = nu
    if "eftw0wa" in outs and hasattr(cp, "EFTCAMB"):
        meta["eft_w0wa"] = list(cp.EFTCAMB.get_effective_w0wa())


def main(cfg):
    sys.path.insert(0, cfg["tree"])
    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(cfg["tree"]).resolve()), camb.__file__
    args = cfg["args"]
    out: dict[str, np.ndarray] = {}
    meta: dict = {"camb_version": camb.__version__, "errors": {}}
    if cfg["task"] == "generic":
        cp = camb.read_ini(cfg["ini"], no_validate=True)
        if args.get("restore_eft") and hasattr(cp, "EFTCAMB"):
            pass
        setup(camb, cp, args, meta)
        outputs(camb, cp, args, out, meta)
    elif cfg["task"] == "setparams":
        params = dict(args["params"])
        if args.get("eftcamb_class"):
            cp = camb.CAMBparams()
            cp.set_classes(dark_energy_model="EFTCAMB")
            cp = camb.set_params(cp=cp, **params)
        else:
            cp = camb.set_params(**params)
        setup(camb, cp, args, meta)
        outputs(camb, cp, args, out, meta)
    else:
        raise ValueError(cfg["task"])
    out["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(cfg["out"], **out)
    print(json.dumps({"ok": True, "errors": meta["errors"]}))


if __name__ == "__main__":
    try:
        main(json.loads(sys.argv[1]))
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        sys.exit(3)
