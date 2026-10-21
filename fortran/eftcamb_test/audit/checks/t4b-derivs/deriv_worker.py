"""Track t4b-derivs worker: evaluate the EFTCAMB timestep cache on finite-difference stencils.

Runs inside a subprocess that imports camb from one frozen tree only.  For a model ini it computes the CAMB results,
then calls ``EFTCAMB.get_eft_functions`` on stencil points a_c * exp(j h), j = -4..4, for several step sizes h in
ln a, around centres a_c log-spaced in [a_lo, a_hi].  Optionally it also calls ``get_scale_evolution`` on a stencil
(for the fields that are only filled inside ``derivs``: Einstein-equation factors, pi, delta_phi coefficients).

Usage (normally via deriv_check.py):
    python deriv_worker.py '<json cfg>'
cfg keys: tree, ini, out, a_lo, a_hi, n_c, hs, evol_k, evol_a_lo, evol_a_hi, evol_n_c, evol_h, python_set
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

OFFS = np.arange(-4, 5)


def main(cfg: dict) -> None:
    sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
    import harness

    harness._exit_with_parent()
    sys.path.insert(0, cfg["tree"])
    import camb

    assert Path(camb.__file__).resolve().is_relative_to(Path(cfg["tree"]).resolve()), camb.__file__
    t0 = time.time()
    cp = camb.read_ini(cfg["ini"], no_validate=True)
    for path, value in cfg.get("python_set", {}).items():
        obj = cp
        parts = path.split(".")
        for p in parts[:-1]:
            obj = getattr(obj, p)
        setattr(obj, parts[-1], value)
    ini_keys = harness.parse_ini(cfg["ini"])
    if cp.EFTCAMB.EFTflag == 0 and int(ini_keys.get("EFTflag", "0")) != 0:
        # pre-rebase trees (B1) do not read EFTCAMB keys in read_ini
        cp.EFTCAMB.initialize_parameters(cp, ini_keys, print_header=False)
    pert_turn_on_input = float(cp.EFTCAMB.EFTCAMB_pert_turn_on)
    # get_background runs SetParams (EFTCAMB init, RGR, stability) and thermo, without the Cl/transfer integration.
    results = camb.get_results(cp) if cfg.get("full_results") else camb.get_background(cp)
    meta = {
        "EFTflag": int(cp.EFTCAMB.EFTflag),
        "designer": bool(cp.EFTCAMB.EFTCAMB_model_is_designer),
        "pert_turn_on_input": pert_turn_on_input,
        "pert_turn_on_results": float(results.Params.EFTCAMB.EFTCAMB_pert_turn_on),
        "back_turn_on": float(results.Params.EFTCAMB.EFTCAMB_back_turn_on),
        "omk": float(cp.omk),
        "num_nu_massive": int(cp.num_nu_massive),
        "time_results": time.time() - t0,
    }
    try:
        meta["model_name"] = cp.EFTCAMB.model_name()
    except Exception as e:  # noqa: BLE001
        meta["model_name"] = repr(e)
    out: dict[str, np.ndarray] = {}
    hs = np.asarray(cfg["hs"], dtype=float)
    xc = np.linspace(np.log(cfg["a_lo"]), np.log(cfg["a_hi"]), int(cfg["n_c"]))
    x = xc[None, None, :] + hs[:, None, None] * OFFS[None, :, None]  # (n_h, 9, n_c)
    a_all = np.exp(x).ravel()
    t1 = time.time()
    # conformal_time (TimeOfzArr) may require monotonic redshifts: evaluate on the sorted (decreasing z) list
    # (and errors out on exact duplicates at z > ~1e10): evaluate on the unique sorted list, then map back
    a_u, inv = np.unique(a_all, return_inverse=True)
    fields, arr_sorted = cp.EFTCAMB.get_eft_functions(results, a_u)
    meta["time_eft_functions"] = time.time() - t1
    meta["fields"] = list(arr_sorted.dtype.names)
    for f in arr_sorted.dtype.names:
        v = np.ascontiguousarray(arr_sorted[f]).ravel()[inv]
        out[f"st_{f}"] = v.reshape(x.shape)
    out["xc"] = xc
    out["hs"] = hs
    # CAMB background at the centres (conformal Hubble from CAMB's dtauda):
    ac = np.exp(xc)
    zc = 1.0 / ac - 1.0
    out["camb_H"] = np.asarray(results.hubble_parameter(zc))  # km/s/Mpc
    out["camb_tau"] = np.asarray(results.conformal_time(zc))
    out["camb_c_km_s"] = np.array(camb.constants.c / 1e3)
    # optional scale evolution on a stencil (single h):
    if cfg.get("evol_k"):
        xe = np.linspace(np.log(cfg["evol_a_lo"]), np.log(cfg["evol_a_hi"]), int(cfg["evol_n_c"]))
        he = float(cfg["evol_h"])
        xev = (xe[None, :] + he * OFFS[:, None]).T.ravel()  # ordered (centre, offset)
        t2 = time.time()
        # get_evolution sorts the times (overlapping stencils interleave): call on the sorted list, then unsort.
        order = np.argsort(xev, kind="stable")
        fields_e, arr_e = cp.EFTCAMB.get_scale_evolution(results, cfg["evol_k"], np.exp(xev[order]))
        meta["time_evolution"] = time.time() - t2
        for f in arr_e.dtype.names:
            v = np.ascontiguousarray(arr_e[f]).reshape(len(cfg["evol_k"]), len(xev))
            u = np.empty_like(v)
            u[:, order] = v
            out[f"ev_{f}"] = u.reshape(len(cfg["evol_k"]), len(xe), len(OFFS))
        a_ev = out["ev_a"].ravel()
        okm = a_ev > 0
        fields_s, arr_s = cp.EFTCAMB.get_eft_functions(results, np.where(okm, a_ev, 0.5))
        for f in arr_s.dtype.names:
            out[f"evst_{f}"] = np.ascontiguousarray(arr_s[f]).reshape(out["ev_a"].shape)
        out["ev_xc"] = xe
        out["ev_h"] = np.array(he)
        out["ev_k"] = np.asarray(cfg["evol_k"], dtype=float)
        out["ev_x_req"] = xev.reshape(len(xe), len(OFFS))
    meta["time_total"] = time.time() - t0
    out["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(cfg["out"], **out)
    print(json.dumps({"ok": True, **{k: v for k, v in meta.items() if k != "fields"}}))


if __name__ == "__main__":
    main(json.loads(sys.argv[1]))
