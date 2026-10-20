"""Verifier for t3-polterddot-curvature-gr: attribution with runtime switches (one diagnostic build of a copy of A0).

DIAGNOSTIC patch (not a proposed fix), selected by env AUDIT_POLTERK (bitmask, read once per process; threads=1):
  bit 1: multiply the octgdot term of EFTCAMB's pigddot by EV%Kf(2)   (as CAMB's denlk2(2) = 3/5 k Kf(2))
  bit 2: multiply the Edot(3) term of EFTCAMB's E2ddot by EV%Kf(2)    (as CAMB's original polterddot formula)
  bit 4: dump k, tau, polter, polterdot, polterddot(EFTCAMB form), 3/50 k octgdot, 1/5 k Edot(3), Kf(2), TC flag
         for k <= AUDIT_POLTERK_KMAX to file AUDIT_POLTERK_FILE (source-time samples only)
AUDIT_POLTERK=0 must reproduce A0 bit for bit.
Usage: python polterK_attrib.py build | run | dump
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

WORK = harness.AUDIT_ROOT / "work" / "verify" / "polterK"
DEST = WORK / "attrib_tree"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"

EDITS = [
    (
        "    real(dl) :: pigddot, E2ddot\n",
        "    real(dl) :: pigddot, E2ddot\n"
        "    integer, save :: audit_pk = -1, audit_pk_unit = -1\n"
        "    real(dl), save :: audit_pk_kmax = 0._dl\n"
        "    character(len=512) :: audit_pk_env\n"
        "    real(dl) :: audit_kfa, audit_kfb, audit_base\n",
    ),
    (
        "            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*octgdot - dopacity*(pig - polter) - &\n"
        "                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot\n"
        "            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*Edot(3)\n"
        "            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot\n",
        "            if (audit_pk < 0) then\n"
        "                call get_environment_variable('AUDIT_POLTERK', audit_pk_env)\n"
        "                audit_pk = 0\n"
        "                if (len_trim(audit_pk_env) > 0) read(audit_pk_env, *) audit_pk\n"
        "                if (iand(audit_pk, 4) /= 0) then\n"
        "                    call get_environment_variable('AUDIT_POLTERK_KMAX', audit_pk_env)\n"
        "                    read(audit_pk_env, *) audit_pk_kmax\n"
        "                    call get_environment_variable('AUDIT_POLTERK_FILE', audit_pk_env)\n"
        "                    open(newunit=audit_pk_unit, file=trim(audit_pk_env), status='replace')\n"
        "                end if\n"
        "            end if\n"
        "            audit_kfa = 1._dl\n"
        "            audit_kfb = 1._dl\n"
        "            if (iand(audit_pk, 1) /= 0) audit_kfa = EV%Kf(2)\n"
        "            if (iand(audit_pk, 2) /= 0) audit_kfb = EV%Kf(2)\n"
        "            pigddot = 2._dl/5._dl*k*qgdot - 3._dl/5._dl*k*audit_kfa*octgdot - dopacity*(pig - polter) - &\n"
        "                opacity*(pigdot - polterdot) + 8._dl/15._dl*k*sigmadot\n"
        "            E2ddot = -dopacity*(E(2) - polter) - opacity*(Edot(2) - polterdot) - k/3._dl*audit_kfb*Edot(3)\n"
        "            polterddot = pigddot/10._dl + 3._dl/5._dl*E2ddot\n"
        "            if (iand(audit_pk, 4) /= 0 .and. k <= audit_pk_kmax) then\n"
        "                audit_base = polterddot + (audit_kfa - 1._dl)*3._dl/50._dl*k*octgdot &\n"
        "                    + (audit_kfb - 1._dl)*k/5._dl*Edot(3)\n"
        "                write(audit_pk_unit, '(9ES25.16E3)') k, tau, polter, polterdot, audit_base, &\n"
        "                    3._dl/50._dl*k*octgdot, k/5._dl*Edot(3), EV%Kf(2), merge(1._dl, 0._dl, EV%TightCoupling)\n"
        "            end if\n",
    ),
]


def build():
    import shutil

    if DEST.exists():
        shutil.rmtree(DEST)
    shutil.copytree(harness.tree_path("A0"), DEST, symlinks=True)
    f = DEST / "fortran" / "equations.f90"
    s = f.read_text()
    for old, new in EDITS:
        assert s.count(old) == 1, old[:60]
        s = s.replace(old, new)
    f.write_text(s)
    print(json.dumps(harness.build_tree(DEST, DEST), indent=1))


OV = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "1500", "do_lensing": "F"}
OMKS = ("0", "-0.05", "0.05", "-0.01", "-0.1", "0.1")
MASKS = (0, 1, 2, 3)


def job(a):
    tree, omk, mask = a
    name = f"{tree if tree in ('A0', 'B0') else 'P'}_m{mask}_omk{omk}"
    return harness.run_model(tree, GR, WORK / "runs" / f"{name}.npz", overrides={**OV, "omk": omk},
                             python_set={"Reion.max_redshift": 50.0}, threads=1,
                             extra_env={"AUDIT_POLTERK": str(mask)})


def stats(a, b):
    a, b = a["cl_unlensed_scalar"], b["cl_unlensed_scalar"]
    r = np.abs(a[2:1500, 0] / b[2:1500, 0] - 1)
    te = np.abs(a[2:1500, 3] - b[2:1500, 3]) / np.sqrt(b[2:1500, 0] * b[2:1500, 1])
    ee = np.abs(a[2:1500, 1] / b[2:1500, 1] - 1)
    return r[:28].max(), 2 + int(np.argmax(r[:28])), r[98:].max(), te.max(), ee.max()


def run():
    jobs = [(t, o, 0) for o in OMKS for t in ("A0", "B0")] + [(str(DEST), o, m) for o in OMKS for m in MASKS]
    with ThreadPoolExecutor(4) as ex:
        for st in ex.map(job, jobs):
            if st["status"] != "ok":
                print("FAILED", st["out"], st["status"], st.get("log_tail"))
    L = lambda n: harness.load(WORK / "runs" / f"{n}.npz")  # noqa: E731
    print(f"{'omk':>6s} {'variant':>22s} {'TT l<30':>10s} {'at':>3s} {'TT l>=100':>10s} {'TE/sqrt':>10s} {'EE':>10s}")
    res = {}
    for o in OMKS:
        b0 = L(f"B0_m0_omk{o}")
        rows = [("A0", L(f"A0_m0_omk{o}"))] + [(f"P mask{m}", L(f"P_m{m}_omk{o}")) for m in MASKS]
        for lab, d in rows:
            s = stats(d, b0)
            res[f"{o}:{lab}"] = s
            print(f"{o:>6s} {lab + ' vs B0':>22s} {s[0]:10.3e} {s[1]:3d} {s[2]:10.3e} {s[3]:10.3e} {s[4]:10.3e}")
        a0 = L(f"A0_m0_omk{o}")["cl_unlensed_scalar"]
        p0 = L(f"P_m0_omk{o}")["cl_unlensed_scalar"]
        print(f"{o:>6s} {'P mask0 == A0 bitwise':>22s} {bool(np.array_equal(a0, p0))}")
    (WORK / "attrib_summary.json").write_text(json.dumps(res, indent=1))


def dump():
    """Numerical consistency: d(polterdot)/dtau along the evolved solution vs the analytic polterddot variants."""
    out = WORK / "dump"
    out.mkdir(parents=True, exist_ok=True)
    for omk in ("-0.05", "0.05"):
        f = out / f"dump_omk{omk}.txt"
        if not f.exists():
            st = harness.run_model(str(DEST), GR, out / f"dump_omk{omk}.npz", overrides={**OV, "omk": omk},
                                   python_set={"Reion.max_redshift": 50.0}, threads=1,
                                   extra_env={"AUDIT_POLTERK": "4", "AUDIT_POLTERK_KMAX": "2e-3",
                                              "AUDIT_POLTERK_FILE": str(f)})
            print(omk, st["status"])
        analyze(f, omk)


def analyze(f, omk):
    d = np.loadtxt(f)
    ks = np.unique(d[:, 0])
    print(f"omk={omk}: {len(d)} rows, {len(ks)} k values <= 2e-3, kmin={ks.min():.4e}")
    X, Y, rows = [], [], []
    for kk in ks:
        s = d[(d[:, 0] == kk) & (d[:, 8] == 0)]
        s = s[np.argsort(s[:, 1])]
        _, iu = np.unique(s[:, 1], return_index=True)
        s = s[iu]
        if len(s) < 50:
            continue
        tau, pdot, base, A, B, kf = s[:, 1], s[:, 3], s[:, 4], s[:, 5], s[:, 6], s[:, 7]
        num = np.gradient(pdot, tau, edge_order=2)
        sel = slice(5, -5)
        scale = np.median(np.abs(num[sel]))
        rel = [np.median(np.abs(num[sel] - c[sel])) / scale
               for c in (base, base + (1 - kf) * A, base + (1 - kf) * B, base + (1 - kf) * (A + B))]
        rows.append((kk, kf[0], *rel))
        if abs(kf[0]) > 1e-12:
            X.append(np.column_stack([(1 - kf[sel]) * A[sel], (1 - kf[sel]) * B[sel]]) / scale)
            Y.append((num[sel] - base[sel]) / scale)
    rows = np.array(rows)
    labs = ("EFTCAMB as is (no Kf)", "Kf on octg only", "Kf on E3 only", "Kf on both (CAMB)")
    for i, lab in enumerate(labs):
        r = rows[:, 2 + i]
        print(f"  {lab:>22s}: median|d(polterdot)/dtau - cand|/median|.| over k: min {r.min():.2e} "
              f"median {np.median(r):.2e} max {r.max():.2e}")
    z = rows[np.abs(rows[:, 1]) < 1e-12]
    for row in z:
        print(f"  Kf2=0 mode k={row[0]:.4e}: " + ", ".join(f"{lab} {v:.2e}" for lab, v in zip(labs, row[2:])))
    X, Y = np.vstack(X), np.concatenate(Y)
    c, *_ = np.linalg.lstsq(X, Y, rcond=None)
    print(f"  LSQ (modes with Kf2/=0) coefficient on (1-Kf2)*3/50 k octgdot = {c[0]:.5f}, "
          f"on (1-Kf2)*1/5 k Edot3 = {c[1]:.5f}   (1 = term needs Kf(2), 0 = EFTCAMB form consistent)")


if __name__ == "__main__":
    {"build": build, "run": run, "dump": dump}[sys.argv[1]]()
